#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 12:31:52 2025

@author: nive
"""

#Import Modules
import pandas as pd

# Load quantification table (feature table) from coverm
# This table contains read counts for CDS
# Choose between:
#   - Inactivation genes: 
#       - unbinned_proteins_inactivation_CDS_ARGs_feature_table.tsv
#       - Blastp_inactivation_protein_CDS_mapping_dataset1.tsv
#   - Efflux genes: 
#       - unbinned_proteins_efflux_CDS_ARGs_feature_table.tsv
#       - Blastp_efflux_protein_CDS_mapping_dataset1.tsv
table1 = pd.read_csv("../../../Data/Dataset1/unbinned_proteins_efflux_CDS_ARGs_feature_table.tsv", sep="\t")  # Quantification table

# Load BLAST mapping table (annotation) 
# Contains CARD protein hits, descriptions, and nucleotide/CDS info
table2 = pd.read_csv("../../../Data/Dataset1/Blastp_efflux_protein_CDS_mapping_dataset1.tsv", sep="\t")  # Annotation table

# Merge feature table with annotation table
# Merge on CDS_ID (from annotation) vs Contig (from feature table)
# 'inner' merge ensures we only keep CDS present in the quantification table
merged = pd.merge(table1, table2, left_on="Contig", right_on="CDS_ID", how="inner")

# Identify read count columns automatically
# Assumes columns ending with "Read Count" contain sample read counts
read_count_cols=[]
for col in merged.columns:
    if col.endswith("Read Count"):
        read_count_cols.append(col)

# Define final columns to keep: annotation + read counts
final_cols = ["CARD_Protein_Accession", "CARD_Protein_Description", "NCBI_Nucleotide_ID","CDS_Name"] + read_count_cols

# Subset merged table to final columns
final = merged[final_cols]

# Save the combined table to TSV under Data/Dataset1 as this will be used as an input file for downstream analysis
# The filename will be either:
#   - inactivation_CDS_read_count_table_dataset1.tsv
#   - efflux_CDS_read_count_table_dataset1.tsv
# depending on which input files were used
final.to_csv("../../../Data/Dataset1/efflux_CDS_read_count_table_dataset1.tsv", sep="\t", index=False)