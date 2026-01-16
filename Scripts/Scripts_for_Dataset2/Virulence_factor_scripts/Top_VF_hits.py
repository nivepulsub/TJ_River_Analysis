#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 16:34:53 2025

@author: nive
"""
# Import modules
import pandas as pd
import re

# Input file containing raw VF search results for each bin
# Each Bin_ID contains multiple Query_ID entries and corresponding Hit_IDs
file_path = "../../../Data/Dataset2/virulence_factors_for_every_bin_dataset2.txt"

records = []      # Stores parsed records for output
record = {}       # Temporary storage for current hit entry
current_query = None

# Parse the raw text file line-by-line
# The logic below extracts: Bin_ID (MAG identifier), Query_ID, Hit_ID, Description, Score, Position

with open(file_path, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        # Detect start of a new MAG bin
        if line.startswith("Bin_ID:"):
            # If we already have a record in progress, save it before starting a new bin
            if record:
                records.append(record)
                record = {}
            record["Bin_ID"] = line.split(":", 1)[1].strip()
        
        # Detect start of a new Query
        elif line.startswith("Query.ID:"):
            # When a new Query appears, save the previous record
            if record.get("Query_ID") and record.get("Hit_ID"):
                records.append(record)
                record = {"Bin_ID": record["Bin_ID"]}
            record["Query_ID"] = line.split(":", 1)[1].split("\t")[0].split(";")[0].strip()
        
        # Extract the first hit only (top hit)
        elif line.startswith("Hit.ID:"):
            # If we already have a hit, skip additional hits for the same Query_ID
            if "Hit_ID" in record:
                continue    # skip additional hits for same query
            parts = line.split("\t")
            record["Hit_ID"] = parts[0].split(":", 1)[1]
            record["Hit_Description"] = parts[1].split(":", 1)[1] if len(parts) > 1 else ""
            record["Subject_Length"] = parts[2].split(":", 1)[1] if len(parts) > 2 else ""
        
        # Extract score-related fields
        elif line.startswith("Score:"):
            fields = line.split("\t")
            for field in fields:
                key, val = field.split(":", 1)
                record[key.strip()] = val.strip()
        
        # Extract alignment start/end positions for query/subject
        elif line.startswith("Query_Start_Position:"):
            fields = line.split("\t")
            for field in fields:
                key, val = field.split(":", 1)
                record[key.strip()] = val.strip()

# Append the last record
if record:
    records.append(record)

# Convert parsed list of dicts into DataFrame
df = pd.DataFrame(records)

# Keep only relevant annotation columns
columns_to_keep = [
    'Bin_ID','Query_ID','Hit_ID','Hit_Description','Subject_Length',
    'Score','Bits','E-value','Identities','Positives','Gaps',
    'Query_Start_Position','Query_End_Position','Subject_Start_Position','Subject_End_Position'
]
df = df[columns_to_keep]

# Keep only the top hit per Bin_ID + Query_ID
df_top_hits = df.drop_duplicates(subset=['Bin_ID','Query_ID'], keep='first')

# Save to TSV
df_top_hits.to_csv("../../../Output/Dataset2/VFs_parsed_top_hits_dataset2.tsv", sep="\t", index=False)

#Load the TSV file
df = pd.read_csv("../../../Output/Dataset2/VFs_parsed_top_hits_dataset2.tsv", sep="\t")


#  Extract VF category from hit descriptions
def extract_vf_category(hit_description):
    matches = re.findall(r"-\s*([A-Za-z0-9\s/:-]+?)\s*\(VFC\d+\)", str(hit_description))
    if matches:
        # take first match
        return matches[0].strip()
    else:
        return "Unclassified"

df["VF_Category"] = df["Hit_Description"].apply(extract_vf_category)

# Normalize text for consistency
df["VF_Category"] = df["VF_Category"].str.strip().str.title()

# Summarize VF counts per MAG (Bin_ID)
vf_summary = df.groupby(["Bin_ID", "VF_Category"]).size().unstack(fill_value=0)
vf_summary["Total_VFs"] = vf_summary.sum(axis=1)

# Save final MAG VF-Count table
vf_summary.to_csv("../../../Output/Dataset2/VF_count_table_per_MAG_dataset2.tsv", sep="\t")


