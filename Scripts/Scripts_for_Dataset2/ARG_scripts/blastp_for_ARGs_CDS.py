# This script was run on a departmental Linux server
# NOTE:
# This script is reused for both efflux and inactivation ARG analyses.
# To switch between efflux and inactivation modes, update:
#   1. query_file
#   2. output_file
#   3. blast_database
#   4. output filenames
# accordingly.

# IMPORT MODULES
import subprocess
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO

# Email for NCBI access
Entrez.email = "nive9756@gmail.com"

# Parameters for BLASTP search against efflux resistance proteins derived from CARD
query_file = "/data/Nive/TJ_dataset2/unbinned_proteins.faa"
output_file = "/data/Nive/TJ_dataset2/blastp_unbinned_proteins_against_efflux_CARD_database.xml"  # Update filename when running inactivation ARG analysis eg., blastp_unbinned_proteins_against_inactivation_CARD_database.xml
blast_executable = "/data/Nive/TJ_dataset2/blast_env/ncbi-blast-2.16.0+/bin/blastp"
blast_database = "/data/Nive/TJ_dataset2/blast_db_antibiotic_efflux_protein/ARG_efflux" # Update the database when running inactivation ARG analysis eg., blast_db_antibiotic_inactivation_protein/ARG_inactivation
evalue = 10                         # E-value threshold
query_coverage_hsp_percentage = 100 # only consider hits covering 100% of the query protein
max_target_seqs = 5                 # Save top 5 hits per query
outfmt = 5                          # XML output format (for Biopython parsing)
num_threads = 36                    # Number of threads for parallel BLAST

# RUN BLAST
# Purpose: Compare unbinned protein sequences to a CARD-derived ARG database using BLASTP
blastp_cmd = [
    blast_executable,
    "-query", query_file,
    "-db", blast_database,
    "-out", output_file,
    "-evalue", str(evalue),
    "-qcov_hsp_perc", str(query_coverage_hsp_percentage),
    "-max_target_seqs", str(max_target_seqs),
    "-outfmt", str(outfmt),
    "-num_threads", str(num_threads)
]

#Run the blastp_cmd
subprocess.run(blastp_cmd, check=True)

# PARSE BLAST XML output
# Writes a text file with query ID, hit IDs, descriptions, alignment scores, percent identities/positives/gaps, alignment positions, and GenBank records.
def parse_blast_output(xml_file, txt_output):
    with open(xml_file) as file: #xml_file: Path to BLAST XML output
        output_parsing = list(NCBIXML.parse(file))
    with open(txt_output, "w") as fin: #txt_output: Path to text file to save parsed results
        for record in output_parsing:
            if record.alignments:
                fin.write("Query.ID:"+record.query+"\tQuery_Length:"+str(record.query_length)+"\n\n")
                for alignment in record.alignments:
                    print("Hit Id:", alignment.hit_id)
                    print("Hit Description:", alignment.hit_def)
                    fin.write("Hit.ID:"+alignment.hit_id+"\tHit_Description:"+alignment.hit_def+"\tSubject_Length:"+str(alignment.length)+"\n\n")
                    for hsp in alignment.hsps:
                        fin.write("Score:"+str(hsp.score)+"\tBits:"+str(hsp.bits)+"\tE-value:"+str(hsp.expect)+"\tIdentities:"+str(round((hsp.identities/hsp.align_length)*100))+"%"+"\tPositives:"+str(round((hsp.positives/hsp.align_length)*100))+"%"+"\tGaps:"+str(round((hsp.gaps/hsp.align_length)*100))+"%"+"\n\n")
                        fin.write("Query_Start_Position:"+str(hsp.query_start)+"\tQuery_End_Position:"+str(hsp.query_end)+"\tSubject_Start_Position:"+str(hsp.sbjct_start)+"\tSubject_End_Position:"+str(hsp.sbjct_end)+"\n\n")
                        fin.write("Query:"+hsp.query+"\n")
                        fin.write("      "+hsp.match+"\n")
                        fin.write("Sbjct:"+hsp.sbjct+"\n\n\n")
                    try:
                        # Fetch GenBank record for hit protein
                        # Sometimes fetching GenBank records is not possible if the accession is not found
                        genbank_id = alignment.hit_def.split('|')[1]
                        print("Fetching:", genbank_id)
                        handle = Entrez.efetch(db="protein", id=genbank_id, rettype="gb", retmode="text")
                        gb_text = handle.read()
                        handle.close()
                        fin.write("----- GenBank Full Record for " + genbank_id + " -----\n")
                        fin.write(gb_text + "\n")
                        fin.write("----- End of GenBank Record -----\n\n")
                    except Exception as e:
                        # We catch exceptions and write an error instead of stopping the script
                        print("Error fetching:", genbank_id, str(e))
                        fin.write("GenBank_Description: Error fetching (" + str(e) + ")\n")

# Call the parse_blast_output function
parse_blast_output(output_file, "Blastp_parsing_results_unbinned_proteins_efflux_ARGs.txt")

# Fetch CDS from protein accession
def fetch_cds_from_protein(protein_acc, fasta_out):
    
    try:
        # Step 1: Link protein accession to nuccore (nucleotide) record
        link_handle = Entrez.elink(dbfrom="protein", db="nuccore", id=protein_acc, linkname="protein_nuccore")
        link_record = Entrez.read(link_handle)
        link_handle.close()

        if not link_record[0]["LinkSetDb"]:
            print(f"No nuccore link found for {protein_acc}")
            return

        nuccore_id = link_record[0]["LinkSetDb"][0]["Link"][0]["Id"]

        # Step 2: Fetch GenBank record to locate CDS region
        nuccore_handle = Entrez.efetch(db="nuccore", id=nuccore_id, rettype="gb", retmode="text")
        gb_record = SeqIO.read(nuccore_handle, "genbank")
        nuccore_handle.close()

        # Step 3: Find the CDS feature that matches the protein accession
        for feature in gb_record.features:
            if feature.type == "CDS" and "protein_id" in feature.qualifiers:
                if protein_acc in feature.qualifiers["protein_id"]:
                    cds_seq = feature.extract(gb_record.seq)
                    cds_record = SeqIO.SeqRecord(
                        cds_seq,
                        id=protein_acc,
                        description=f"CDS corresponding to protein {protein_acc} from {gb_record.id}"
                    )
                    # Step 4: Write CDS sequence to THE FASTA file
                    with open(fasta_out, "a") as fout:
                        SeqIO.write(cds_record, fout, "fasta")
                    print(f"CDS sequence written for {protein_acc}")
                    return  # Stop after first match

        print(f"No CDS found for {protein_acc} in {nuccore_id}")

    except Exception as e:
        print(f"Error fetching CDS for {protein_acc}: {e}")

# Extract top protein hits
def extract_top_protein_hits(xml_file, protein_out, nucleotide_out):

    with open(xml_file) as file:
        blast_records = list(NCBIXML.parse(file))

    with open(protein_out, "w") as fout_protein:
        for record in blast_records:
            if not record.alignments:
                print(f"No hits found for {record.query}")
                continue

            top_alignment = record.alignments[0]

            try:
                hit_acc = top_alignment.hit_def.split('|')[1]
            except IndexError:
                print(f"Could not parse accession from: {top_alignment.hit_def}")
                continue

            hit_desc = top_alignment.hit_def
            print(f"Fetching top protein hit for {record.query}: {hit_acc} - {hit_desc}")

            try:
                # Get the protein sequence
                handle = Entrez.efetch(db="protein", id=hit_acc, rettype="fasta", retmode="text")
                seq_record = SeqIO.read(handle, "fasta")
                handle.close()

                seq_record.id = f"{record.query}|{hit_acc}"
                seq_record.description = hit_desc
                SeqIO.write(seq_record, fout_protein, "fasta")

                print(f"Wrote protein sequence for {hit_acc}")

                # Fetch corresponding CDS sequence
                fetch_cds_from_protein(hit_acc, nucleotide_out)

            except Exception as e:
                print(f"Error fetching sequences for {hit_acc}: {e}")

extract_top_protein_hits(
    "blastp_unbinned_proteins_against_efflux_CARD_database.xml", # full parsed BLAST report # Update filename when running inactivation ARG analysis eg., blastp_unbinned_proteins_against_inactivation_CARD_database.xml
    "Blastp_top_efflux_hits_from_unbinned_proteins.faa", # protein sequences of top hits # Update filename when running inactivation ARG analysis eg., Blastp_top_inactivation_hits_from_unbinned_proteins.faa
    "Blastp_top_efflux_hits_CDS_from_unbinned_proteins.fna" # CDS sequences corresponding to top protein hits  # Update filename when running inactivation ARG analysis eg., Blastp_top_inactivation_hits_CDS_from_unbinned_proteins.fna
)
