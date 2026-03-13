library(ggplot2)
library(Seurat)
library(dplyr)
library(data.table)
library(future)
library(DoubletFinder)  # Ensure DoubletFinder is installed
options(future.globals.maxSize = 50000 * 1024^2)
# Set up parallel computing
plan(multicore, workers = 80)

# Find all h5 files ending with "_cb_filtered_seurat.h5" in current directory
h5_files <- list.files(pattern = "*_cb_filtered_seurat\\.h5$", full.names = TRUE)

# If no files found with the specified pattern, try other h5 file patterns
if (length(h5_files) == 0) {
  h5_files <- list.files(pattern = "*.h5$", full.names = TRUE)
  cat("No files ending with '_cb_filtered_seurat.h5' found, but found the following h5 files:\n")
  print(h5_files)
}

if (length(h5_files) == 0) {
  stop("No h5 files found in current directory")
}

cat("Found", length(h5_files), "h5 files\n")

# Process each file in a for loop
for (h5_file in h5_files) {
  # Extract sample name: remove "_cb_filtered_seurat.h5" suffix
  sample_name <- gsub("_cb_filtered_seurat\\.h5$", "", basename(h5_file))
  cat("\n===================================\n")
  cat("Processing sample:", sample_name, "\n")
  cat("File path:", h5_file, "\n")
  cat("===================================\n")
  
  tryCatch({
    # 1. Read h5 file
    cat("Reading h5 file...\n")
    data <- Read10X_h5(h5_file)
    cat("Reading complete. Data dimensions:", dim(data)[1], "genes x", dim(data)[2], "cells\n")
    
    # 2. Create Seurat object
    cat("Creating Seurat object...\n")
    seurat_obj <- CreateSeuratObject(
      counts = data, 
      min.cells = 3, 
      min.features = 0,
      project = sample_name
    )
    
    # 3. Add sample information
    seurat_obj$Sample <- sample_name
    
    # 4. Quality filtering
    cat("Filtering low-quality cells (nCount_RNA >= 500 & nFeature_RNA >= 1000)...\n")
    seurat_obj <- subset(seurat_obj, subset = nCount_RNA >= 500 & nFeature_RNA >= 1000)
    cat("Number of cells after filtering:", ncol(seurat_obj), "\n")
    
    # 5. Normalization
    cat("Normalizing data...\n")
    seurat_obj <- NormalizeData(
      object = seurat_obj, 
      normalization.method = "LogNormalize", 
      scale.factor = 1e4
    )
    
    # 6. Find variable features
    cat("Finding variable features...\n")
    seurat_obj <- FindVariableFeatures(object = seurat_obj)
    
    # 7. Scale data
    cat("Scaling data...\n")
    seurat_obj <- ScaleData(
      object = seurat_obj, 
      features = rownames(x = seurat_obj), 
      vars.to.regress = c("nCount_RNA")
    )
    
    # 8. PCA dimensionality reduction
    cat("Performing PCA dimensionality reduction...\n")
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    
    # 9. UMAP dimensionality reduction
    cat("Performing UMAP dimensionality reduction...\n")
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE, reduction = "pca")
    
    # 10. Clustering
    cat("Performing clustering...\n")
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
    
    # 11. Doublet detection
    cat("Detecting doublets...\n")
    sweep.res.list <- paramSweep(seurat_obj, PCs = 1:10, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    
    annotations <- seurat_obj@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(0.075 * ncol(seurat_obj@assays$RNA))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    
    seurat_obj <- doubletFinder(
      seurat_obj, 
      PCs = 1:20, 
      pN = 0.25, 
      pK = mpK, 
      nExp = nExp_poi, 
      reuse.pANN = NULL, 
      sct = FALSE
    )
    
    # 12. Save doublet detection results
    colnames(seurat_obj@meta.data)[length(seurat_obj@meta.data)] <- "DoubleScore"
    doublet <- table(seurat_obj$DoubleScore)  
    table_name <- paste0(sample_name, "_doublet_cell_number.txt")
    write.table(doublet, table_name, quote = FALSE, sep = "\t", row.names = FALSE)
    cat("Doublet detection results saved to:", table_name, "\n")
    
    # 13. Filter out doublets
    seurat_obj <- subset(seurat_obj, cells = WhichCells(seurat_obj, expression = `DoubleScore` != "Doublet"))
    cat("Number of cells after removing doublets:", ncol(seurat_obj), "\n")
    
    # 14. Save results
    output_file <- paste0(sample_name, "_singlet.rds")
    saveRDS(seurat_obj, output_file)
    cat("Processing complete! Results saved as:", output_file, "\n")
    
    # 15. Clean up memory
    rm(data, seurat_obj)
    gc()
    
  }, error = function(e) {
    cat("Error processing file:", h5_file, "\n")
    cat("Error message:", e$message, "\n")
    cat("Skipping this file, continuing to next...\n")
  })
}

cat("\n===================================\n")
cat("All files processed!\n")
cat("===================================\n")
