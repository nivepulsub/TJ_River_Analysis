# This script was run on a departmental Linux server
# Import modules
import subprocess
from Bio.Blast import NCBIXML
from Bio import Entrez
import os
import glob

# Parameters for blastp search against VFDB (Virulence Factor Database) 
query_folder = "/mnt/HDD/niveps/binning_output_faa_files/"  # Folder with MAG protein files (.faa)
output_file = "/mnt/HDD/niveps/blastp_against_virulence_factors_against_VFDB_dataset1.xml"
output_folder = "/mnt/HDD/niveps/binned_genomes_folder_to_find_VF_dataset1"
blast_executable = "/mnt/HDD/niveps/blast_env/ncbi-blast-2.16.0+/bin/blastp"
blast_database = "/mnt/HDD/niveps/blast_db_virulence_factor/virulence_factor"
evalue = 10                               # E-value cutoff
query_coverage_hsp_percentage = 100       # Only consider hits covering 100% of the query protein
max_target_seqs = 5                       # Save top 5 hits per query
outfmt = 5                                # XML output format (for Biopython parsing)
num_threads = 36                          # Number of threads for parallel BLAST


# Parse the .xml output file
def parse_blast_output(xml_file,output_txt,bin_id):
   # Parse XML using Biopython
   with open(xml_file) as file:
      output_parsing = list(NCBIXML.parse(file))
   
   # Append parsed results to the output text file
   with open(output_txt, "a") as fin:
      Entrez.email = "nive9756@gmail.com"

# Loop through each BLAST records to get the details of query and subject
      for record in output_parsing:
         if record.alignments:
            fin.write("Bin_ID:"+bin_id+"\n\n")
            fin.write("Query.ID:"+record.query+"\tQuery_Length:"+str(record.query_length)+"\n\n")
            for alignment in record.alignments:
               print("Hit Id:",alignment.hit_id)
               print("Hit Description:",alignment.hit_def)
               fin.write("Hit.ID:"+alignment.hit_id+"\tHit_Description:"+alignment.hit_def+"\tSubject_Length:"+str(alignment.length)+"\n\n")
               for hsp in alignment.hsps:
                  fin.write("Score:"+str(hsp.score)+"\tBits:"+str(hsp.bits)+"\tE-value:"+str(hsp.expect)+"\tIdentities:"+str(round((hsp.identities/hsp.align_length)*100))+"%"+"\tPositives:"+str(round((hsp.positives/hsp.align_length)*100))+"%"+"\tGaps:"+str(round((hsp.gaps/hsp.align_length)*100))+"%"+"\n\n")
                  fin.write("Query_Start_Position:"+str(hsp.query_start)+"\tQuery_End_Position:"+str(hsp.query_end)+"\tSubject_Start_Position:"+str(hsp.sbjct_start)+"\tSubject_End_Position:"+str(hsp.sbjct_end)+"\n\n")
                  fin.write("Query:"+hsp.query+"\n")
                  fin.write("      "+hsp.match+"\n")
                  fin.write("Sbjct:"+hsp.sbjct+"\n\n\n")
               # Try to fetch GenBank annotation for the matched protein
               try:
                  genbank_id = alignment.hit_def.split('|')[1].split(')')[0]
                  print("Fetching:", genbank_id)

                  handle = Entrez.efetch(db="protein", id=genbank_id, rettype="gb", retmode="text")
        
                  # Write entire GenBank record as plain text
                  gb_text = handle.read()
                  handle.close()

                  fin.write("----- GenBank Full Record for " + genbank_id + " -----\n")
                  fin.write(gb_text + "\n")
                  fin.write("----- End of GenBank Record -----\n\n")

               except Exception as e:
                  print("Error fetching:", genbank_id, str(e))
                  fin.write("GenBank_Description: Error fetching (" + str(e) + ")\n")


# Define a function called run_blastp for multiple runs
def run_blastp(query_file, output_file, blast_database):
# Command for blastp
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

# Run the BLAST command
   subprocess.run(blastp_cmd,check=True)

# Process each .faa file
for query_file in glob.glob(os.path.join(query_folder, "*.faa")):
   base_name = os.path.splitext(os.path.basename(query_file))[0]
   output_file = os.path.join(output_folder, f"{base_name}_vfdb.xml")
# Call the run_blastp function
   run_blastp(query_file,output_file,blast_database)
# Call the parse_blast_output function
   parse_blast_output(output_file,"virulence_factors_for_every_bin_dataset1.txt",base_name)

