create_peak_matrix_IDX <- function(seurat_rds, outpath) {

  # Load required packages
  library(data.table)
  library(dplyr)
  library(Seurat)

  # Step 1: Load Seurat object and extract peak data
  message("Step 1: Loading Seurat object and extracting peaks...")
  seurat_obj <- readRDS(seurat_rds)
  peaks <- GetAssayData(seurat_obj, assay = assay_name, slot = slot_name)
  peak_names <- rownames(peaks)

  IDX <- data.frame(peaks = paste0(seqnames,"-",ranges), idx = as.integer(1:length(peak_names)))
  saveRDS(IDX, paste0(outpath,"NHPABC_allpeaks_idx.rds"))

  return(IDX)

}
