# **Analysis of epigenetic erosin in NHPABC**

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
