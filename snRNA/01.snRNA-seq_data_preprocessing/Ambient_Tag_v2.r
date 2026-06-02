#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(Seurat)
library(tidyr)
library(tidyverse)
library(DropletQC)
library(patchwork)
library(magrittr)
library(parallel)
library(GeneOverlap)
library(utils)
library(ggpubr)
source('/path/to/custom/functions.R')  # Load custom functions

# Set global parameters
options(repr.plot.width = 12, repr.plot.height = 12)
timesExcess = 2
res = 0.5
npcs = 30
ncores = 8
ambCutoffCoef = 1.5
logfc = 0.25

# Get mb parameter from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript process_mb.R <mb>")
}
mb <- args[1]

# Construct file paths
target_dir <- file.path("/path/to/input", mb)
results_dir <- file.path("/path/to/output", mb)
annotation_path <- "/path/to/annotation/genomic.gtf"

# Set version and output directory names
version <- "v2"
txt_dir <- paste0(mb, "_count.csv")
sub_dirs <- c('M4_Str_1', 'M4_Str_2', 'M5_Str_1', 'M5_Str_2')

# Main processing loop
for (sub_dir in sub_dirs) {
  cat("Processing:", sub_dir, "\n")
  
  # Set up directories
  raw_data_dir <- file.path(target_dir, sub_dir)
  res_sub_dir <- file.path(results_dir, sub_dir)
  dir.create(res_sub_dir, recursive = TRUE)
  figure_output_dir <- file.path(res_sub_dir, "Figure")
  dir.create(figure_output_dir, recursive = TRUE)
  
  # STEP 1: Read raw data
  raw_data <- Read10X(data.dir = file.path(raw_data_dir, "FilterMatrix"), gene.column = 1)
  sampObj <- CreateSeuratObject(counts = raw_data)
  num_rows1 <- ncol(sampObj)
  
  # STEP 2-3: Filter and cluster
  sampObj <- subset(sampObj, subset = nCount_RNA > 0)
  finalBarcs <- read.table(gzfile(file.path(raw_data_dir, "FilterMatrix/barcodes.tsv.gz")), row.names = 1)
  
  # Label cells as empty/non-empty
  sampObj$is_empty <- ifelse(colnames(sampObj) %in% rownames(finalBarcs), 'Non_empty', 'Empty')
  
  # Soft filtering
  if (prop.table(table(sampObj$is_empty))['Empty'] < 0.5) {
    cat("More than half of the cell barcodes were not filtered. Keeping all barcodes.\n")
    softFiltObj <- sampObj
  } else {
    rawMeta <- sampObj@meta.data
    rawMetaFinalCount <- sum(rawMeta$is_empty == 'Non_empty')
    keepMax <- min(rawMetaFinalCount * timesExcess, nrow(rawMeta))
    keepL <- rawMeta[order(rawMeta$nCount_RNA, decreasing = TRUE)[1:keepMax], ] %>% rownames
    softFiltObj <- subset(sampObj, cells = keepL)
  }
  
  # Run clustering
  cat("Running clustering with excess barcodes\n")
  softFiltObj <- NormalizeData(softFiltObj, verbose = FALSE)
  softFiltObj <- FindVariableFeatures(softFiltObj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  softFiltObj <- ScaleData(softFiltObj, features = rownames(softFiltObj), verbose = FALSE)
  softFiltObj <- RunPCA(softFiltObj, verbose = FALSE)
  softFiltObj <- FindNeighbors(softFiltObj, dims = 1:npcs, verbose = FALSE, reduction = 'pca')
  softFiltObj <- FindClusters(softFiltObj, verbose = FALSE, resolution = res)
  
  # Identify ambient clusters
  filtObjMeta <- softFiltObj[[]]
  clDF <- data.frame(clusters = names(table(filtObjMeta$seurat_clusters)), 
                     size = as.numeric(table(filtObjMeta$seurat_clusters)))
  clDF$emptySize <- as.numeric(table(filtObjMeta[filtObjMeta$is_empty == 'Empty', 'seurat_clusters']))
  clDF$ratio <- clDF$emptySize / clDF$size
  
  cutoff <- prop.table(table(sampObj$is_empty))['Empty'] * ambCutoffCoef
  clDF_Ambient <- clDF[clDF$ratio > cutoff | clDF$size == max(clDF$size), ]
  
  filtObjMeta$is_ambient_cluster <- 'Non_Ambient'
  filtObjMeta[filtObjMeta$seurat_clusters %in% clDF_Ambient$clusters, 'is_ambient_cluster'] <- 'Ambient'
  softFiltObj@meta.data <- filtObjMeta
  sampObjAmb <- softFiltObj
  sampObjAmb <- RunUMAP(sampObjAmb, dims = 1:30, verbose = FALSE)
  
  # Save figures before cellbender
  pdf(file.path(figure_output_dir, "beforeCB_ambCluster.pdf"), width = 12)
  p <- DimPlot(sampObjAmb, group.by = 'is_ambient_cluster', label = TRUE, raster = TRUE, label.size = 7) +
    theme(text = element_text(size = 20, face = 'bold')) + NoLegend() + ggtitle('')
  print(p)
  dev.off()
  
  pdf(file.path(figure_output_dir, "beforeCB_DotPlot.pdf"), width = 12)
  p <- DotPlot(sampObjAmb, features = c('SNAP25', 'SLC17A7', 'GAD1', 'GAD2', 'PLP1', 'MBP', 'SLC1A3', 'GFAP', 'CSF1R', 'A2M', 'PDGFRA')) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  dev.off()
  
  # STEP 5-6: Process cellbender output and annotate
  file_path <- file.path(raw_data_dir, "cb_filtered_seurat.h5")
  hmat <- Read10X_h5(filename = file_path, use.names = TRUE)
  seurat_obj <- CreateSeuratObject(counts = hmat)
  
  # Add metadata
  sample <- strsplit(sub_dir, "_")[[1]][1]
  region <- strsplit(sub_dir, "_")[[1]][2]
  library_id <- strsplit(sub_dir, "_")[[1]][3]
  seurat_obj <- AddMetaData(seurat_obj, 
                            metadata = data.frame(Sample = rep(sample, ncol(seurat_obj)),
                                                  Library = rep(library_id, ncol(seurat_obj)),
                                                  Region = rep(region, ncol(seurat_obj)),
                                                  batch = rep(version, ncol(seurat_obj))))
  
  # Standard Seurat processing
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.2, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:15, verbose = FALSE)
  
  # Transfer ambient cluster labels
  raw_meta <- sampObjAmb@meta.data
  new_meta <- seurat_obj@meta.data
  common_indices <- intersect(rownames(raw_meta), rownames(new_meta))
  new_meta[common_indices, "is_ambient_cluster"] <- raw_meta[common_indices, "is_ambient_cluster"]
  seurat_obj <- AddMetaData(seurat_obj, metadata = new_meta$is_ambient_cluster, col.name = "is_ambient_cluster")
  
  # Save figures after cellbender
  pdf(file.path(figure_output_dir, "afterCB_ambCluster.pdf"), width = 12)
  p <- DimPlot(seurat_obj, group.by = 'is_ambient_cluster', label = TRUE, raster = TRUE, label.size = 7) +
    theme(text = element_text(size = 20, face = 'bold')) + NoLegend() + ggtitle('')
  print(p)
  dev.off()
  
  pdf(file.path(figure_output_dir, "afterCB_DotPlot.pdf"), width = 12)
  p <- DotPlot(seurat_obj, features = c('SNAP25', 'SLC17A7', 'GAD1', 'GAD2', 'PLP1', 'MBP', 'SLC1A3', 'GFAP', 'CSF1R', 'A2M', 'PDGFRA')) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  dev.off()
  
  # Save RDS
  saveRDS(seurat_obj, file.path(res_sub_dir, "cb_ambTag.RDS"))
  
  # Save summary statistics
  csv_path <- file.path(results_dir, txt_dir)
  meta1 <- sampObjAmb@meta.data
  ambient_counts1 <- table(meta1$is_ambient_cluster)
  meta2 <- seurat_obj@meta.data
  num_rows2 <- nrow(meta2)
  ambient_counts2 <- table(meta2$is_ambient_cluster)
  
  if (file.exists(csv_path)) {
    existing_data <- read.csv(csv_path, row.names = 1)
  } else {
    existing_data <- data.frame(
      bfCB_cells = integer(),
      bfCB_Ambient = integer(),
      bfCB_Non_Ambient = integer(),
      afCB_cells = integer(),
      afCB_Ambient = integer(),
      afCB_Non_Ambient = integer(),
      stringsAsFactors = FALSE
    )
  }
  
  new_sample <- data.frame(
    bfCB_cells = num_rows1,
    bfCB_Ambient = as.integer(ambient_counts1["Ambient"]),
    bfCB_Non_Ambient = as.integer(ambient_counts1["Non_Ambient"]),
    afCB_cells = num_rows2,
    afCB_Ambient = as.integer(ambient_counts2["Ambient"]),
    afCB_Non_Ambient = as.integer(ambient_counts2["Non_Ambient"]),
    stringsAsFactors = FALSE
  )
  rownames(new_sample) <- sub_dir
  existing_data <- rbind(existing_data, new_sample)
  write.csv(existing_data, file = csv_path, row.names = TRUE)
  
  cat("Process completed:", sub_dir, "\n")
}
