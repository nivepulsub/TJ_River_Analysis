# This script was run on a departmental Linux server
# NOTE:
# This script is reused for both efflux and inactivation ARG analyses.
# To switch between efflux and inactivation modes, update:
#   1. BLAST XML input file
#   2. Output protein FASTA filename
#   3. Output CDS FASTA filename
#   4. Output mapping TSV filename
#   5. CARD database used during the BLASTP step

# Import Modules
from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO

def extract_top_protein_hits(xml_file, protein_out, cds_out, mapping_out):
    Entrez.email = "nive9756@gmail.com"

    # Parse BLAST XML file
    with open(xml_file) as file:
        blast_records = list(NCBIXML.parse(file))

    # Output files
    with open(protein_out, "w") as fout_protein, open(mapping_out, "w") as fmap:
        fmap.write("Query_ID\tCARD_Protein_Accession\tCARD_Protein_Description\tNCBI_Nucleotide_ID\tCDS_Name\tLocus_Tag\tCDS_ID\n")

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

            nucleotide_id = "NA"
            cds_name = "NA"
            locus_tag = "NA"
            cds_id = "NA"

            try:
                # Fetch protein FASTA
                handle = Entrez.efetch(db="protein", id=hit_acc, rettype="fasta", retmode="text")
                seq_record = SeqIO.read(handle, "fasta")
                handle.close()

                seq_record.id = f"{record.query}|{hit_acc}"
                seq_record.description = hit_desc
                SeqIO.write(seq_record, fout_protein, "fasta")
                print(f"Wrote protein sequence for {hit_acc}")

                # Find linked nucleotide record
                link_handle = Entrez.elink(dbfrom="protein", db="nuccore", id=hit_acc, linkname="protein_nuccore")
                link_record = Entrez.read(link_handle)
                link_handle.close()

                if link_record[0]["LinkSetDb"]:
                    nucleotide_id = link_record[0]["LinkSetDb"][0]["Link"][0]["Id"]

                    # Fetch GenBank file
                    gb_handle = Entrez.efetch(db="nuccore", id=nucleotide_id, rettype="gb", retmode="text")
                    gb_record = SeqIO.read(gb_handle, "genbank")
                    gb_handle.close()
                    nucleotide_id = gb_record.id

                    # Look for CDS feature matching protein_id
                    for feature in gb_record.features:
                        if feature.type == "CDS" and "protein_id" in feature.qualifiers:
                            if hit_acc in feature.qualifiers["protein_id"]:
                                cds_id = feature.qualifiers["protein_id"][0]
                                cds_name = feature.qualifiers.get("gene", ["NA"])[0]
                                locus_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]

                                # Extract CDS sequence
                                cds_seq = feature.extract(gb_record.seq)
                                cds_record = SeqIO.SeqRecord(
                                    cds_seq,
                                    id=f"{record.query}|{hit_acc}|{cds_id}",
                                    description=f"CDS for protein {hit_acc} ({cds_name}) from {gb_record.id}"
                                )
                                with open(cds_out, "a") as fout_cds:
                                    SeqIO.write(cds_record, fout_cds, "fasta")

                                print(f"Wrote CDS for {hit_acc}: {cds_name} ({cds_id})")
                                break
                    else:
                        print(f"No CDS found for {hit_acc} in {nucleotide_id}")

            except Exception as e:
                print(f"Error fetching info for {hit_acc}: {e}")

            # Write mapping summary
            fmap.write(f"{record.query}\t{hit_acc}\t{hit_desc}\t{nucleotide_id}\t{cds_name}\t{locus_tag}\t{cds_id}\n")

# NOTE: Filenames below correspond to inactivation ARG analysis.
# Update filenames accordingly when running efflux ARG analysis.
# (e.g., replace "inactivation" with "efflux" in input and output filenames).

if __name__ == "__main__":
    extract_top_protein_hits(
        "blastp_unbinned_proteins_against_inactivation_CARD_database.xml",
        "Blastp_top_CDS_inactivation_hits_from_unbinned_proteins_run1.faa",
        "Blastp_top_CDS_inactivation_hits_from_unbinned_proteins_run1.fna",
        "Blastp_inactivation_protein_CDS_mapping_dataset1.tsv"
    )
