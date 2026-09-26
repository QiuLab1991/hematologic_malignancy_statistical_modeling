This repository contains the core computational code supporting the manuscript “Prospective Plasma Proteogenomics Identifies Pre-diagnostic Protein Signatures Associated with Hematologic Malignancies.” 

**Script	Main analysis**

**01_proteomics_cox.R**	Prospective plasma-proteomic association analyses using Cox proportional-hazards models, including the full cohort, region-based internal replication, random 60:40 robustness analysis, and Benjamini-Hochberg FDR correction.

**02_pseudotime_mfuzz.R**	Group-level pseudo-temporal analysis of pre-diagnostic protein profiles. Protein abundance is summarized across time-to-diagnosis bins, standardized, and clustered using fuzzy c-means clustering implemented in Mfuzz.

**03_genetic_mr.R**	Genetic analyses using harmonized protein instruments and hematologic-malignancy GWAS summary statistics. Wald-ratio and fixed-effect inverse-variance-weighted estimates are generated as appropriate, with heterogeneity, MR-Egger intercept, and FDR analyses.

**04_scrna_virtual_knockout.R**	Exploratory FL single-cell analyses, including candidate-gene expression across FL and reference lymph-node samples, B-cell subset analyses, and virtual knockout of FAS, ICAM1, and LTBR in FL tumor B cells using scTenifoldKnk.


**Software**

The statistical analyses were performed using R version 4.4.2. Major R packages used in the scripts include:
- data.table
- survival
- dplyr
- Biobase
- Mfuzz
- TwoSampleMR
- ggplot2
- Seurat
- scTenifoldKnk


**Data**

The repository does not redistribute individual-level or controlled-access datasets. UK Biobank data are available to approved researchers through the UK Biobank access procedures. Publicly available pQTL, GWAS, eQTL, and single-cell resources and their accession information are described in the manuscript and its Supplementary Data.


**Citation**

Wei C, Li Y, Luo Y, Li H, Zhou X, Yang L, Wei Q, Niu T, Chen Z, Qiu S. Prospective Plasma Proteogenomics Identifies Pre-diagnostic Protein Signatures Associated with Hematologic Malignancies. 2026. Code repository available at: https://github.com/QiuLab1991/hematologic_malignancy_statistical_modeling.
