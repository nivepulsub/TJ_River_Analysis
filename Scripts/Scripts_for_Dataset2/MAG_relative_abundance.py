#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 14:08:59 2025

@author: nive
"""
# Import modules
import pandas as pd

# Load MAG relative abundance table
# genome coverage table with relative abundance per MAG per sample
df = pd.read_csv("../../Data/Dataset2/genome_coverage_table_dataset2.tsv", sep="\t")

# Remove the row where Genome == 'unmapped'
df = df[df['Genome'].str.lower() != 'unmapped']

# Keep only bacterial MAGs
df_gtdb = pd.read_csv("../../Data/Dataset2/gtdbtk.bac120.summary.tsv", sep="\t")
bacterial_bins = df_gtdb['user_genome'].tolist()
df = df[df['Genome'].isin(bacterial_bins)]

# Keep only relative abundance columns
relative_abundance_cols = [col for col in df.columns if "Relative Abundance (%)" in col]
print(relative_abundance_cols)

# Compute total relative abundance per MAG
# Fill missing values with 0 and sum across all samples
df['TotalAbundance'] = df[relative_abundance_cols].fillna(0).sum(axis=1)

#Sort MAGs by total abundance (descending)
df_sorted = df.sort_values(by='TotalAbundance', ascending=False)

# Keep only Genome, Relative Abundance columns, and TotalAbundance
columns_to_keep = ['Genome'] + relative_abundance_cols + ['TotalAbundance']
final_table = df_sorted[columns_to_keep]

# Save to CSV
final_table.to_csv("../../Output/Dataset2/MAGs_relative_abundance_dataset2.csv", sep=",", index=False)

# Load the CSV file
df = pd.read_csv("../../Output/Dataset2/MAGs_relative_abundance_dataset2.csv", sep=",")

# Create feature table: sampleIDs = rows, MAGs(Genomes) = columns
# Subset to Genome + relative abundance cols
feature_table = df[['Genome'] + relative_abundance_cols]

# Transpose so that sampleIDs are rows and MAGs are columns 
# Set Genome as index and transpose so that samples become rows
feature_table = feature_table.set_index('Genome').T
print(feature_table)

# Clean up row names (remove "Relative Abundance (%)" to keep sample IDs clean)
feature_table.index = feature_table.index.str.replace(" Relative Abundance (%)", "", regex=False)
feature_table.index = feature_table.index.str.replace(".sorted", "", regex=False)

# Save feature table
feature_table.to_csv("../../Output/Dataset2/MAGs_feature_table_dataset2.csv",index=True, index_label="SampleID")

# Load the feature table and Taxonomy file
df1 = pd.read_csv("../../Output/Dataset2/MAGs_feature_table_dataset2.csv", sep=",") 
df2 = pd.read_csv("../../Data/Dataset2/gtdbtk.bac120.summary.tsv", sep="\t")

# Extract Bin_ID for MAG_ID creation
df2["Bin_ID"] = df2["user_genome"].str.split(".").str[1]  # "binning_output.001" → "001"

# Create MAG_ID (Choose S1 for Study1 or S2 for Study2 depending on the study) 
df2["MAG_ID"] = "S2-" + df2["Bin_ID"]

# Create a dictionary for mapping
user_to_tax = df2.set_index("user_genome")["classification"].to_dict()
user_to_mag = df2.set_index("user_genome")["MAG_ID"].to_dict()


# Create map for rank labels
rank_suffix = {
    "s__": "(S)",
    "g__": "(G)",
    "f__": "(F)",
    "o__": "(O)",
    "c__": "(C)",
    "p__": "(P)",
    "d__": "(D)"
}

def get_lowest_rank(tax_string):
    tax_string = str(tax_string).strip()

    # Handle unclassified
    if tax_string.lower().startswith("unclassified bacteria"):
        return "Unclassified Bacteria (D)"
    if tax_string.lower().startswith("unclassified"):
        return "Unclassified (U)"

    ranks = tax_string.split(";")

    # species -> genus -> family -> ...
    for tax in reversed(ranks):
        for prefix, tag in rank_suffix.items():
            if tax.startswith(prefix):
                name = tax.replace(prefix, "").strip()
                if name != "":
                    return f"{name} {tag}"

    return tax_string  


# Now rename each column: turn full taxonomy into shortened taxonomy + MAG_ID
new_columns = []
for col in df1.columns:
    if col == "SampleID":
        new_columns.append(col)
        continue

    user_genome = col  # e.g., binning_output.001

    # Force an error if user_genome is not found
    try:
        full_tax = user_to_tax[user_genome]
        mag_id = user_to_mag[user_genome]
    except KeyError:
        raise KeyError(f"Column '{user_genome}' not found in GTDB file. "
                       f"Check if your df1 column names match df2['user_genome'].")

    short_tax = get_lowest_rank(full_tax)
    # Always include MAG_ID
    if mag_id != "":
        new_label = f"{short_tax} {mag_id}"
    else:
        new_label = short_tax
    new_columns.append(new_label)

# Assign new column names
df1.columns = new_columns

# Save final table under Data/Dataset2 as this will serve as an input file for NMDS analysis
df1.to_csv("../../Data/Dataset2/MAGs_feature_table_lowest_taxonomy_dataset2.csv", index=False)