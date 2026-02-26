#!/bin/bash

# Directory containing input MAG FASTA files
INPUT_DIR="binning_output_fasta_files"

# Directory to store predicted protein files
OUTPUT_DIR="binning_output_faa_files"

# Create output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

# Loop through each FASTA file and run Prodigal
for file in ${INPUT_DIR}/*.fasta; do
    base=$(basename "$file" .fasta)

    # Run Prodigal in metagenomic mode
    prodigal -i "$file" \
             -a "${OUTPUT_DIR}/${base}.faa" \
             -p meta
done