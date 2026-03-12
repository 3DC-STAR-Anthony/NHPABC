# **snATAC-seq Data Processing Pipeline**
**Overview**: This pipeline processes single-nucleus ATAC-seq (snATAC-seq) data from raw sequencing reads to annotated cell clusters, with quality control and integration with matched snRNA-seq markers.

## Step 1: Read Alignment & Preprocessing
- Tool: PISA (https://github.com/shiquan/PISA)
- Reference: T2T-MFA8, 2025, V1.1 genome
- Output: Fragment files for each library (compatible with ArchR)

## Step 2: Cell Quality Control (QC)
- Filtering Criteria:
  - TSS enrichment score < 4
  - Fragment count < 3,000
- Doublet Removal:
  - Doublet score calculated with addDoubletScores(ArchR)
  - Filtered with filterDoublets(filterRatio = 2)
```r
#Generate_arrow file
source("~/Create_arrow_files.R")
arrow_files <- create_arrow_files(
  input_dir = input_dir,
  genome_annotation_path = "/path/to/T2TMF8_genomeAnnotation.rds",
  gene_annotation_path = "/path/to/T2TMF8_geneAnnotation.rds",
  minTSS = 4,
  minFrags = 3000,
  threads = 40
)
```
## Step 3: Dimensionality Reduction & Clustering
- Feature Creation: 500-bp genomic tiles
- Dimensionality Reduction: Iterative Latent Semantic Indexing (LSI) via addIterativeLSI(ArchR)
- Clustering: Seurat algorithm applied to LSI dimensions (resolution = 0.8)

## Step 4: Cluster Annotation
- Method: Annotated using marker genes from matched snRNA-seq dataset
- Integration: Consistent annotation between snATAC-seq and snRNA-seq modalities

