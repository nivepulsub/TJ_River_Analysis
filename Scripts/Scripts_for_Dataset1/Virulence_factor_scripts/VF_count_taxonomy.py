#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 16:39:31 2025

@author: nive
"""
# Import modules
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load VF count table and taxonomy
VF_count_table = pd.read_csv("../../../Output/Dataset1/VF_count_table_per_MAG_dataset1.tsv", sep="\t")
Taxonomy_output = pd.read_csv("../../../Data/Dataset1/gtdbtk.bac120.summary.tsv", sep="\t")

# Select only columns we need from Taxonomy_output
taxonomy = Taxonomy_output[["user_genome", "classification"]].rename(
    columns={"user_genome": "Bin_ID", "classification": "Taxonomy"}
)

# Merge VF counts with taxonomy
merged = pd.merge(VF_count_table, taxonomy, on="Bin_ID", how="left")

# Convert to string and strip whitespace
merged["Taxonomy"] = merged["Taxonomy"].astype(str).str.strip()

# Remove empty strings and NaNs
merged = merged[merged["Taxonomy"].notna()]
merged = merged[merged["Taxonomy"] != ""]
merged = merged[merged["Taxonomy"].str.lower() != "nan"]


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

def get_lowest_rank(x):
    x = x.strip()

    # For Unclassified Bacteria
    if x.lower().startswith("unclassified bacteria"):
        return "Unclassified Bacteria (D)"
    if x.lower().startswith("unclassified"):
        return "Unclassified (U)"
    
    ranks = x.split(";")
    
    # Reverse order: species -> genus -> family -> order -> class -> phylum -> domain
    for tax in reversed(ranks):
        for prefix, tag in rank_suffix.items():
            if tax.startswith(prefix):
                name = tax.replace(prefix, "").strip()
                if name != "":
                    return f"{name}{tag}"
    
# Keep the original full taxonomy for reference
merged["Full_Taxonomy"] = merged["Taxonomy"].copy()

# Apply the function
merged["Taxonomy"] = merged["Taxonomy"].apply(get_lowest_rank)

# Create MAG_ID for Study 1: S1-001, S1-002,... (or) For study 2: S2-001, S2-002,... 
# Bin_ID format is 'binning_output.001','binning_output.002',...... for both the studies
# Replace S1 or S2 depending on which study dataset we are working with
merged["MAG_ID"] = "S1-" + merged["Bin_ID"].str.split('.').str[1]

# Combine Taxonomy and MAG_ID for heatmap row labels
merged["Taxonomy_Label"] = merged["Taxonomy"] + " " + merged["MAG_ID"]

# Reorder VF columns and put Others last
vf_cols = [c for c in merged.columns if c not in ["Bin_ID", "Taxonomy", "Total_VFs"]]

if "Others" in vf_cols:
    vf_cols.remove("Others")
    vf_cols.append("Others")

# VF acronyms
vf_acronyms = {
    "Adherence": "Adherence",
    "Antimicrobial Activity/Competitive Advantage": "Antimicrobial",
    "Biofilm": "Biofilm",
    "Effector Delivery System": "EDS",
    "Exoenzyme": "Exoenzyme",
    "Exotoxin": "Exotoxin",
    "Immune Modulation": "Immune Mod.",
    "Invasion": "Invasion",
    "Motility": "Motility",
    "Nutritional/Metabolic Factor": "Nutri/Metab",
    "Post-Translational Modification": "PTM",
    "Regulation": "Regulation",
    "Stress/Survival": "Stress/Survival",
    "Others": "Others"
}

# Rename the columns
merged = merged.rename(columns=vf_acronyms)

# VF order
vf_order = [
    "Adherence", "Antimicrobial", "Biofilm", "EDS", "Exoenzyme",
    "Exotoxin", "Immune Mod.", "Invasion", "Motility", "Nutri/Metab",
    "PTM", "Regulation", "Stress/Survival", "Others"
]

# Use the acronym names in the same order
vf_cols = [vf for vf in vf_order if vf in merged.columns]

# Sort by Total VFs
merged["Total_VFs"] = merged[vf_cols].sum(axis=1)
merged = merged.sort_values("Total_VFs", ascending=False)

# Save the table
merged.to_csv("../../../Output/Dataset1/VF_count_with_taxonomy_dataset1.tsv", sep="\t", index=False)

# Heatmap for all MAGs
vf_data_all = merged[vf_cols]
vf_data_all.index = merged["Taxonomy_Label"]

plt.figure(figsize=(5, 55))
sns.heatmap(vf_data_all, cmap="Reds")

plt.xticks(rotation=45, ha='right')
plt.title("Virulence Factor Category Abundance (Dataset 1)", fontsize=12)

plt.savefig("../../../Output/Dataset1/VF_counts_heatmap_dataset1_all.pdf",
            bbox_inches='tight', dpi=300)
plt.show()

# Heatmap for Top 20 MAGs
top20 = merged.head(20)

vf_data_20 = top20[vf_cols]
vf_data_20.index = top20["Taxonomy_Label"]

plt.figure(figsize=(5,4))
sns.heatmap(vf_data_20, cmap="Reds")

plt.xticks(rotation=45, ha='right')
plt.title("Virulence Factor Category Abundance (Dataset 1)", fontsize=12)

plt.savefig("../../../Output/Dataset1/VF_counts_heatmap_dataset1_top20.pdf",
            bbox_inches='tight', dpi=300)
plt.show()