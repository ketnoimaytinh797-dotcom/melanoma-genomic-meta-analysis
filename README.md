# Genomic alterations in cutaneous melanoma – analysis code

This repository contains the R analysis script and study-level datasets used for the meta-analysis of genomic alterations in cutaneous melanoma.

## Files

`meta_mutation.R`  
R script used for the statistical analyses reported in the manuscript.

`Table S11 Data from Studies Included in the Meta-analysis.xlsx`  
Study-level dataset used as the primary input for the analyses.

`Table S12 List of Studies and Reasons for Exclusion.xlsx`  
List of excluded studies and reasons for exclusion.

`Asian Melanoma.xlsx`  
Asian subset dataset used to estimate pooled mutation prevalence in Asian populations.

## Overview of the analysis

The analyses were conducted in R using study-level mutation data extracted from eligible studies.

The script performs:
- random-effects meta-analysis of mutation prevalence
- subgroup analyses by gene analysis methodology
- pooled mutation prevalence estimates for primary and metastatic melanoma
- pooled mutation prevalence analyses in Asian populations

## Input data

The primary dataset used in the analysis is:

`Table S11 Data from Studies Included in the Meta-analysis.xlsx`

Subgroup analyses by sequencing methodology are based on the `Methodology` column in this dataset, while study identifiers correspond to the `TT` column.

## Running the analysis

Place the input Excel files in the same directory as `meta_mutation.R`, then run the script in R or RStudio:

```r
source("meta_mutation.R")
```

Output tables and summary files will be written to the `results/` directory.

## Code availability

The analysis code and supporting datasets are shared to improve transparency and reproducibility of the analyses reported in the manuscript.
