# Motif Enrichment Analysis of NHPABC

- Motif Database: [`JASPAR2022`](https://doi.org/10.1093/nar/gkab1113).
- TF Filtering: Excluded TFs expressed in <5% of cells (based on matched snRNA-seq clusters).
- Enrichment Calculation: Computed motif enrichment in target peak sets versus background using ArchR's computeEnrichment.

## Step 1: Preparations of Input File
```r
source("~/Create_motif_peak_matrix.R")
motif_matrix <- create_motif_peak_matrix(
  seurat_rds = "~/PFC_Astrocyte_PeakMatrix.rds",
  genome_annotation_rds = "~/T2TMF8_genomeAnnotation.rds",
  species = "Homo sapiens",
  collection = "CORE",
  p_cutoff = 5e-05,
  motif_width = 7
)
```
```r
source("")
IDX <- create_peak_matrix_IDX(
  seurat_rds = "~/PFC_Astrocyte_PeakMatrix.rds",
  outpath = "~"
)
```
## Step 2:Enrichment Calculation
```r
source("~/Compute_dare_motif_enrichment.R")
# Load data
IDX <- IDX
motifMat <- motif_matrix
dare <- fread("~/NHPABC_dare_final.csv", data.table = FALSE) #pDARE, sDARE-Early et al.

# Run enrichment analysis
results <- compute_dare_motif_enrichment(
  dare = dare,
  motifMat = motifMat,
  IDX = IDX,
  outpath = "/path/to/output/directory"
)
# Access results
up_results <- results$up
down_results <- results$dow
```
## Step 3：TF Filtering
This step was based on the expression results of the TF genes in the corresponding snRNA-seq, which is not elaborated here.
