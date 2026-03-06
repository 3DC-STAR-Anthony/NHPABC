# Identification of candidate cis-regulatory elements(cCREs) in NHPABC
This project performs progressive/stage-specific/longevity-associated differentially accessible regulatory elements (pDAREs/sDAREs/longDAREs) analysis on single-nuclei ATAC-seq data. 
# pDAREs
## Input
- An ArchR PeakMatrix object in `.rds` format named: `<Subtype>.rds` (e.g., `Ast1_PFC_PeakMatrix.rds`)
- Use the function of ArchR [`getMatrixFromProject`](https://www.archrproject.com/reference/getMatrixFromProject.html)
- The object must contain:
  - peakmatrix (peaks X Individuals)
  - Metadata columns including:
    - Annotion (e.g., `Subtype`)
    - Individuals (e.g., `Y1,Y2,..,VO6`)
## Step 1:


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
