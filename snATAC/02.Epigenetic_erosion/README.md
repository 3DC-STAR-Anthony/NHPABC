# **Analysis of epigenetic erosin in NHPABC**
Prior to analysis, cell subtypes were filtered to ensure robust representation.For each cell subtype, individuals with fewer than 10 nuclei were excluded, and only cell subtypes with at least 100 nuclei and more than two individuals per age group were retained.
## TSSEnrichment
```r
source("~/Calculate_epigenetic_erosion_lm.R")
tss_results <- analyze_age_qc_correlation(
  proj = readRDS('/path/to/project.rds'),
  metric = "TSS",
  output_path = "/path/to/TSS_lm_result.csv"
)
```
## FRIP
```r
source("~/Calculate_epigenetic_erosion_lm.R")
frip_results <- analyze_age_qc_correlation(
  proj = readRDS('/path/to/project.rds'),
  metric = "FRIP",
  output_path = "/path/to/FRIP_lm_result.csv"
)
```
