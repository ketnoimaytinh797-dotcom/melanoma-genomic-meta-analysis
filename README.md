# Genomic alterations in cutaneous melanoma: analysis code

This repository contains the R script and study-level data files used for the meta-analysis of genomic alterations in cutaneous melanoma.

## Repository contents

- `meta_mutation.R`  
  Main analysis script. The script reads the study-level workbook directly from `Table S11 Data from Studies Included in the Meta-analysis.xlsx`, standardizes the `Gene mutation` sheet, and runs the pooled prevalence analyses.

- `Table S11 Data from Studies Included in the Meta-analysis.xlsx`  
  Primary study-level dataset used for the meta-analysis.

- `Table S12 List of Studies and Reasons for Exclusion.xlsx`  
  Excluded studies and reasons for exclusion. This file is included for completeness of the systematic review workflow and is not used in the pooled statistical models.

- `Asian Melanoma.xlsx` *(optional)*  
  Asia-specific study-level dataset used for pooled prevalence estimates in the Asian cohort. If this file is not present, the Asia-specific analysis block is skipped automatically.

## Analysis overview

The script performs the following steps:

1. reads the `Gene mutation` sheet from the study-level workbook;
2. reconstructs consistent variable names for study-level metadata, primary-tumor gene columns, and metastatic-tumor gene columns;
3. standardizes the gene analysis method using the `Methodology` field;
4. performs random-effects meta-analyses of mutation prevalence using the `meta` package (`sm = "PFT"`, inverse-variance pooling, DerSimonian–Laird estimator);
5. runs the Table 2 method-based subgroup analyses for the key genes (`BRAF`, `NRAS`, `TP53`, `TERT`, `NF1`);
6. optionally runs the Asia-specific pooled prevalence analyses for primary and metastatic tumors when `Asian Melanoma.xlsx` is available; and
7. writes cleaned analysis inputs, model objects, and plain-text summaries to the `results/` directory.

For the method-based subgroup analyses, study methods are assigned from the `Methodology` column, and study identifiers are reported from the `TT` column.

## Software requirements

The script uses the following R packages:

- `readxl`
- `dplyr`
- `stringr`
- `purrr`
- `readr`
- `meta`

## How to run

1. Place the Excel workbooks in the repository root:
   - `Table S11 Data from Studies Included in the Meta-analysis.xlsx`
   - `Table S12 List of Studies and Reasons for Exclusion.xlsx`
   - `Asian Melanoma.xlsx` *(optional)*
2. Open R or RStudio in the repository directory.
3. Install the required packages if needed.
4. Run:

```r
source("meta_mutation.R")
```

No `setwd()` call is required; the script uses relative paths.

## Output files

Running the script creates a `results/` directory containing:

- cleaned analysis-ready CSV files;
- method-mapping and study-mapping tables;
- `.rds` files with the fitted meta-analysis models; and
- plain-text summaries of the pooled prevalence models.

## Notes

- The primary manuscript dataset is `Table S11 Data from Studies Included in the Meta-analysis.xlsx`.
- The Table 2 method subgroup analyses are defined from the `Methodology` field rather than the free-text `Method` field.
- The script is written so that the main analysis can be reproduced directly from the study-level workbook without intermediate input files.
