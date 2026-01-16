#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 12:41:48 2025

@author: nive
"""

# Import Modules
import pandas as pd
import numpy as np

# Load inactivation ARG read count table
df_1 = pd.read_csv("../../../Data/Dataset2/inactivation_CDS_read_count_table_dataset2.tsv", sep="\t")

# Remove duplicate rows based on CARD_Protein_Accession
# Ensures unique CDS entries for downstream analysis
df_unique_1 = df_1.drop_duplicates(subset=['CARD_Protein_Accession'])

# Load efflux ARG read count table
df_2 = pd.read_csv("../../../Data/Dataset2/efflux_CDS_read_count_table_dataset2.tsv", sep="\t")

# Remove duplicate rows based on CARD_Protein_Accession
df_unique_2 = df_2.drop_duplicates(subset=['CARD_Protein_Accession'])

# Merge inactivation and efflux CDS tables
# Stack them row-wise to create a combined unbinned ARG table
merged = pd.concat([df_unique_1, df_unique_2], ignore_index=True)
# Save to new TSV
merged.to_csv("../../../Output/Dataset2/unbinned_proteins_inactivation_efflux_CDS_ARGs_merged_dataset2.tsv", sep="\t", index=False)

# Summarize Top 10 ARGs by total read count
df = merged.copy()

# Extract gene name from CARD_Protein_Description
# The format is gb|Accession|Gene [Organism]
df['ARG'] = df['CARD_Protein_Description'].str.split('|').str[2].str.split(' ').str[0]

# Keep only read count columns and the category identifier
read_count_cols = [col for col in df.columns if "Read Count" in col]
df["ARG Category"] = np.where(df["CARD_Protein_Accession"].isin(df_unique_1["CARD_Protein_Accession"]),"Inactivation", "Efflux")

# Calculate total read count for each ARG (sum across all samples)
df["Total Read Count"] = df[read_count_cols].sum(axis=1)

# Get Top 10 ARGs by total read count across all samples
top10_ARGs = df.nlargest(10, "Total Read Count")[["ARG","ARG Category","Total Read Count"]]
top10_ARGs["Total Read Count"] = top10_ARGs["Total Read Count"].round(1)

# Save Top 10 ARGs table
top10_ARGs.to_csv("../../../Output/Dataset2/Top10_ARGs_by_Total_Read_Count.tsv", sep="\t", index=False)

# Load the merged table
df = pd.read_csv("../../../Output/Dataset2/unbinned_proteins_inactivation_efflux_CDS_ARGs_merged_dataset2.tsv", sep="\t")
# Keep only the Read Count columns
read_count_cols = [col for col in df.columns if "Read Count" in col]
print(read_count_cols)

# Keep only the CARD_Protein_Accession and SampleID
columns_to_keep = ['CARD_Protein_Accession'] + read_count_cols 
feature_table = df[columns_to_keep]

# Transpose so that sampleIDs are rows and ARGs are columns
feature_table = feature_table.set_index('CARD_Protein_Accession').T
print(feature_table)

# Clean sample names (remove "Read Count" or ".sorted")
feature_table.index = feature_table.index.str.replace("Read Count", "", regex=False)
feature_table.index = feature_table.index.str.replace(".sorted", "", regex=False)

# Remove rows where ALL values are zero
feature_table = feature_table[(feature_table != 0).any(axis=1)]

# Map CARD_Protein_Accession to ARG names
card_info = df.copy()

# Extract only the ARG name from CARD_Protein_Description
card_info["ARG_name"] = (
    card_info["CARD_Protein_Description"]
    .str.split("|").str[2]         # get the part after the second '|'
    .str.split(" \[").str[0]       # split at ' [' and take first part
    .str.strip()                   # remove any leading/trailing spaces
)

# Create mapping: CARD_Protein_Accession -> ARG_name
arg_mapping = card_info.set_index("CARD_Protein_Accession")["ARG_name"].to_dict()

# Rename columns in feature_table
feature_table = feature_table.rename(columns=arg_mapping)

# Saved the feature table under Data/Dataset2 as this will be used as an input file for NMDS analysis
feature_table.to_csv("../../../Data/Dataset2/ARGs_CDS_feature_table_dataset2.csv", index=True, index_label="SampleID")



