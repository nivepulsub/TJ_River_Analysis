#!/bin/bash

# Read file pairs from unaligned_reads_list2
while read file1 file2; do
    # Run fastq_pair for each pair
    fastq_pair "$file1" "$file2"
done < unaligned_reads_list2
