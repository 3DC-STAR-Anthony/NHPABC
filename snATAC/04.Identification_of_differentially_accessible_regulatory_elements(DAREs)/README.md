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
cor <- Process_subtype_cpm_cor(
  input_rds = "~/Ast1_PFC_cCRE.rds",
  output_dir = "~/output",
  mc_cores = 140
)
```
## Step 2: Calculate the coefficient across age
```r
source("~/Process_coefficient.R")
coef <- process_coefficient(
  input_rds = "~/Ast1_PFC_cCRE.rds",
  output_dir = "~/output"
)
```
## Step 3: Filter by corrlation and coefficient
```r
source("~/Filter_correlation_with_coefficient.R")
filter_correlation_with_coefficient(
  subtype_name = "Ast1_PFC",
  cor_results = cor,
  coeff_results = coef,
  output_dir = "~/output"
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
Identification of longevity-associated DAREs (longDAREs)：We used [`getMarkerFeatures`](https://www.archrproject.com/reference/getMatrixFromProject.html?search-input=getMarkerFeatures) function in ArchR to compare the nuclei in the exceptionally old group with the rest of the nuclei for each cell type. cCREs with a Pval < 0.0002 and an absolute log2[fold change] > 1.1 were considered as longDAREs.
## Step 1: Get MarkerFeatures
```r
Rscript Run_getMarkerFeatures.R
```
## Step 2:Filter
```r
scMarkers <- readRDS("Ast1_PFC.LongDAREs_raw.rds")
markerList = getMarkers_new(scMarkers, cutOff = "FDR <= 1")
markerList_filter <- markerList[markerList$Log2FC > 1.1 & markerList$Pval <= 0.0002, ]
```
## Step 3:Add p2g and consvervation information
