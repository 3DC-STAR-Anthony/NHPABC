# **snATAC-seq Data Processing Pipeline**
Overview: This pipeline processes single-nucleus RNA-seq (snRNA-seq) data from raw sequencing reads to annotated cell clusters, with quality control.

## Step 1: Data Processing
- Alignment: Raw reads were aligned to the T2T-MFA8 (v1.1) reference genome using a custom workflow
- Read Counting: Both exonic and intronic reads were included to capture nuclear pre-mRNA
- Doublet Removal: DoubletFinder (v2.0.3) with default settings
- Ambient RNA Removal: CellBender
- Quality Control:
  - Retained nuclei with >500 detected genes and >500 UMIs
  - Excluded nuclei with >1% ribosomal gene content

## Step 2: Downstream Analysis
- Normalization: SCTransform-based method
- Dimensionality Reduction:
  - PCA on top 3,000 highly variable genes
  - Harmony for batch effect correction
- Visualization: UMAP with 30 PCs, 10 neighbors
- Clustering: Manual annotation using canonical cell markers
- Subclustering: Same methodology applied for detailed cell type analysis

# Run processing
```r

```
