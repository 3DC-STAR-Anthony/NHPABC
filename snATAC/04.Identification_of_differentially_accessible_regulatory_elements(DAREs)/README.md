# Identification of candidate cis-regulatory elements(cCREs) in NHPABC
This project performs progressive/stage-specific/longevity-associated differentially accessible regulatory elements (pDAREs/sDAREs/longDAREs) analysis on single-nuclei ATAC-seq data. 
# pDAREs
## Input
- An mean_CPM matrix object in [`_cCRE.rds`](https://github.com/3DC-STAR-Anthony/NHPABC/blob/main/snATAC/03.Identification_of_candidate_cis-regulatory_elements(cCREs)/calculate_cCRE.R) format named: `<Subtype>_cCRE_row.rds` (e.g., `Ast1_PFC_cCRE.rds`)
- The object must contain:
  - mean_CPM matrix (peaks X Individuals)
## Step 1: Calculate the corrlation across age
```r
source("~/Process_subtype_cpm_cor.R")
Process_subtype_cpm_cor(
  input_rds = "~/Ast1_PFC_cCRE.rds",
  output_dir = "~/output",
  mc_cores = 140
)
```
## Step 2:
```r
source("")
result <- process_coefficient_from_rds(
  input_rds = "/path/to/your/subtype_matrix.rds",
  output_dir = "/path/to/output/directory"
)
```
# sDAREs
Use the normalization.method function in Seurat to generate a Peak × Cell CPM (Counts Per Million) matrix for each cell type.
```r
pm <- readRDS("Ast1_PFC_PeakMatrix.rds")
rds <- CreateSeuratObject(counts = pm, assay = "peaks",meta = meta)
rds <- NormalizeData(object = rds, normalization.method = "RC", scale.factor = 1e6) 
saveRDS(rds, paste0("Ast1_PFC_peakmatrix_seurat.rds"))
```
# longDAREs
- Mean CPM > 4 in at least 4 monkey samples;
- Mean CPM > 0 in at least 12 monkey samples.
```r
source("~/calculate_cCRE.R")
cCRE <- calculate_cCRE(
  rds_path = "./Ast1_PFC_peakmatrix_seurat.rds",
  output_dir = "./cCRE_result/"
)
```
