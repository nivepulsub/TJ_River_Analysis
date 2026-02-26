Metagenomic Analysis of Dataset1 and Dataset2:

This project involves datasets from two different studies: Allsing et al., 2022 (ENA Project Accession Number PRJEB57859) and Shahar et al., 2024 (NCBI BioProject Number PRJNA1423341). The input files are stored in the Data/ directory, with separate subfolders for Dataset1 and Dataset2.

All analyses were performed using R, Python, and Jupyter Notebooks, which are organized under the Scripts/ directory, with separate subfolders for Dataset1 and Dataset2, respectively. To facilitate correlation network analysis, the libs/ folder contains the required library, pybootnet.py, which is used to identify associations between antibiotic resistance genes (ARGs) and metagenome-assembled genomes (MAGs).

All outputs are directed to the Output/ directory by default. If an output is required as an input for downstream analyses, it is redirected back to the Data/ directory to maintain proper workflow continuity. Within the Output/ directory, separate subfolders exist for Dataset1 and Dataset2 to keep results organized.

For clarity, the Scripts/Scripts_for_Dataset1/ folder contains two subfolders:
ARG_scripts/: contains scripts for antibiotic resistance gene analysis.
Virulence_factor_scripts/: contains scripts for virulence factor analysis.

Within ARG_scripts/, the following scripts are included:
1. blastp_for_ARGs_CDS.py
2. protein_CDS_mapping.py
3. merging_final_feature_table.py
4. Duplicates_removal_and_merging.py
Some scripts, such as blastp_for_ARGs_CDS.py and protein_CDS_mapping.py, were run on a departmental Linux server. Certain steps may take a long time to execute, and the scripts are configured to use the maximum number of threads to speed up processing. The recommended order of execution is as listed above, starting with blastp_for_ARGs_CDS.py.

Within Virulence_factor_scripts/, the following scripts are included:
1. blastp_for_VFs.py
2. Top_VF_hits.py
3. VF_count_taxonomy.py
Similarly, scripts such as blastp_for_VFs.py were run on a departmental Linux server and may require a significant amount of time to complete. They should be executed in the order listed to maintain workflow continuity.

Bioinformatics Workflow Documentation:
All preprocessing and genome-resolved metagenomic analyses — including quality control, assembly, binning, MAG quality assessment, taxonomic classification, functional annotation, read alignment, and abundance estimation — are documented in detail under the Workflow_Documentation/ directory.

This documentation includes:
1. Software versions
2. Commands executed
3. Parameter explanations
4. Input and output descriptions
5. Workflow numbering corresponding to the workflow diagram (Workflow_overview.pdf)

This separation ensures clarity between:
1. Downstream statistical analyses (Scripts/)
2. Core bioinformatics pipeline steps (Workflow_Documentation/)