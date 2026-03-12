# **Analysis of cell type abundance across aging in NHPABC**
To ensure robust detection, we retained cell types only if they are represented across at least two age groups and ≥2 individuals (each with ≥10 cells)
## Calculate proportions
```r
source("~/Calculate_proportions.R")
meta <- read.csv("~/meta_data.csv")
df <- calculate_proportions(meta)
```
## Identify progressive linear alterations


## Identify stage-specific alterations

