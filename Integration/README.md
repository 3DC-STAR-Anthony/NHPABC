# **Integration of snRNA-seq and snATAC-seq data**

# Integration of global snRNA-seq and snATAC-seq data
For the integration of global snRNA-seq and snATAC-seq datasets, each was randomly downsampled to 300,000 cells.
For specific information, please refer to the document Inregrate_all_30w.R.

# Integration of snRNA-seq and snATAC-seq data in each region
```r
library(ArchR)
library(Seurat)
library(Signac)
library(future)
library(future.apply)
library(parallel)
library(dplyr)
library(RColorBrewer)
library(harmony)
source("~/Integrate_RNA_ATAC.R")
ATAC <- readRDS("~/Save-ArchR-Project.rds")
ATAC_filtered <- ATACflow(ATAC)
RNA <- readRDS("~/PFC_final_annotation_harmony.rds")
rna_atac_integrated <- IntegrateFlow(RNA,ATAC_filtered)
```
