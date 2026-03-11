# Motif Enrichment Analysis of NHPABC

- Motif Database: [`JASPAR2022`](https://doi.org/10.1093/nar/gkab1113).
- TF Filtering: Excluded TFs expressed in <5% of cells (based on matched snRNA-seq clusters).
- Peak Annotation: Annotated peaks with known motifs using ArchR's addMotifAnnotations(default settings).
- Background Selection: Generated GC-matched background peaks for each peak using ArchR's getBgdPeak(default settings).
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
