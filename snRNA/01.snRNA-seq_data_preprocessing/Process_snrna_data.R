library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(future)
library(Matrix)
library(cowplot)
library(pheatmap)
library(Signac)
library(harmony)

process_snrna_data <- function(rds1_path, 
                               rds2_path, 
                               brain_region,
                               output_path,
                               variable_features_n = 3000,
                               npcs = 30,
                               resolution = 0.8,
                               min_pct = 0.1,
                               logfc_threshold = 0.1) {
  
  # Set parallel processing options
  plan(multisession, workers = 1)
  options(future.globals.maxSize = 50000 * 1024^2)
  
  # Load and merge datasets
  message("Loading and merging datasets...")
  rds1 <- readRDS(rds1_path)
  rds2 <- readRDS(rds2_path)
  
  # Add batch information
  rds2$batch <- 'V2'
  
  # Merge datasets
  rds <- merge(rds1, rds2)
  
  # Quality control filtering
  message("Performing quality control filtering...")
  filter <- subset(rds, cells = row.names(rds@meta.data[
    rds@meta.data$DoubleScore != 'Doublet' & 
    rds@meta.data$is_ambient_cluster != 'Ambient', ]))
  
  # Join layers
  filter <- JoinLayers(filter)
  
  # Save raw count data
  raw_count_file <- paste0(output_path, "/", brain_region, "_raw_count.rds")
  message("Saving raw count data: ", raw_count_file)
  saveRDS(filter, raw_count_file)
  
  # Normalization and SCTransform
  message("Normalizing data and performing SCTransform...")
  filter <- NormalizeData(filter, normalization.method = "LogNormalize", scale.factor = 10000)
  filter <- SCTransform(filter, vst.flavor = "v2", verbose = TRUE, 
                        variable.features.n = variable_features_n)
  
  # Filter variable features
  filter@assays$SCT@var.features <- filter@assays$SCT@var.features[
    grep(invert = TRUE, 'LOC1', filter@assays$SCT@var.features)]
  
  # Dimensionality reduction
  message("Performing dimensionality reduction...")
  filter <- RunPCA(filter, npcs = npcs, verbose = TRUE)
  filter <- RunHarmony(filter, group.by.vars = "batch", assay.use = "SCT")
  
  # UMAP and clustering
  message("Running UMAP and clustering...")
  filter <- RunUMAP(filter, reduction = "harmony", dims = 1:npcs, verbose = TRUE)
  filter <- FindNeighbors(filter, reduction = "harmony", dims = 1:npcs, verbose = TRUE)
  filter <- FindClusters(filter, resolution = resolution, verbose = TRUE)
  
  # Find marker genes
  message("Finding marker genes...")
  rds.markers_1 <- FindAllMarkers(object = filter, assay = "SCT", only.pos = TRUE, 
                                   min.pct = min_pct, logfc.threshold = logfc_threshold)
  
  # Save marker genes
  marker_file <- paste0(output_path, "/", brain_region, "_snRNA_findmarker_gene.csv")
  message("Saving marker genes: ", marker_file)
  write.csv(rds.markers_1, marker_file, row.names = FALSE)
  
  # Save final annotated object
  final_file <- paste0(output_path, "/", brain_region, "_final_annotation_harmony.rds")
  message("Saving final annotated object: ", final_file)
  saveRDS(filter, final_file)
  
  # Print summary statistics
  message("\nProcessing complete!")
  message("Brain region: ", brain_region)
  message("Number of cells after filtering: ", ncol(filter))
  message("Number of clusters identified: ", length(unique(filter$seurat_clusters)))
  message("Number of marker genes found: ", nrow(rds.markers_1))
  
  # Create output list
  result <- list(
    seurat_obj = filter,
    markers = rds.markers_1,
  )
  
  return(result)
}
