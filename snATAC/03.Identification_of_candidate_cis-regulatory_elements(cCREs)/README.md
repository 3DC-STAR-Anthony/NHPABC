# Identification of candidate cis-regulatory elements(cCREs) in NHPABC
This project performs cCRE analysis on single-nuclei ATAC-seq data. Cell types containing more than 10 nuclei per individual animal and over 100 nuclei per age group were selected for analysis. cCREs was performed through the following steps:
## Input
- An ArchR PeakMatrix object in `.rds` format named: `<Subtype>.rds` (e.g., `Ast1_PFC_PeakMatrix.rds`)
- Use the function of ArchR [`getMatrixFromProject`](https://www.archrproject.com/reference/getMatrixFromProject.html)
- The object must contain:
  - peakmatrix (peaks X Individuals)
  - Metadata columns including:
    - Annotion (e.g., `Subtype`)
    - Individuals (e.g., `Y1,Y2,..,VO6`)
## Step 1: Convert the peak matrix of each cell type into a Seurat object
Use the normalization.method function in Seurat to generate a Peak × Cell CPM (Counts Per Million) matrix for each cell type.
```r
pm <- readRDS("Ast1_PFC_PeakMatrix.rds")
rds <- CreateSeuratObject(counts = pm, assay = "peaks",meta = meta)
rds <- NormalizeData(object = rds, normalization.method = "RC", scale.factor = 1e6) 
saveRDS(rds, paste0("Ast1_PFC_peakmatrix_seurat.rds"))
```
## Step 2: Further derive a Peak × Individuals mean-CPM matrix, and filter peaks that satisfy
- Mean CPM > 4 in at least 4 monkey samples;
- Mean CPM > 0 in at least 12 monkey samples.
```r
source("~/calculate_cCRE.R")
cCRE <- calculate_cCRE(
  rds_path = "./Ast1_PFC_peakmatrix_seurat.rds",
  output_dir = "./cCRE_result/"
)
```
# Generation of peak-to-gene links
In brief, co-accessible regions were identified for all cCREs in each cell subtype using addCoAccessibility function in ArchR with the following parameters: aggregation k = 10, window size = 500 kb, distance constraint = 250 kb.

# Identification of conserved cCREs
We assessed evolutionary conservation of regulatory elements between monkeys and humans:
- Genomic Conversion: Monkey cCREs were mapped to human genome (hg38) using liftOver
- Orthology Detection: Successfully converted cCREs considered orthologous
- Conservation Validation:
    - Compared against [`human brain snATAC-seq`](https://www.science.org/doi/10.1126/science.adf7044) and [`hippocampus bulk ATAC-seq`](https://www.nature.com/articles/s41586-019-1917-5) data
    - Used bedtools intersection analysis
    - Overlapping cCREs defined as cell-type conserved
## Step 1: Processing human brain and monkey snATAC-seq cCRE data





