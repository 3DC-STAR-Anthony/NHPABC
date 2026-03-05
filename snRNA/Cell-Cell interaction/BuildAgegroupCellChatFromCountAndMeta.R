#library(CellChat)
#library(Matrix)

# Core function to build CellChat objects for all age groups
BuildAgegroupCellChatFromCountAndMeta <- function(base_path, prefix, save_path) {
  # Define age groups and create output directory
  age_groups <- c("Old", "Exceptionally old", "Middle age", "Young")
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
  # Process each age group iteratively
  lapply(age_groups, function(group) {
    # Generate sanitized filenames (replace spaces with underscores)
    sanitized <- gsub(" ", "_", group)
    count_file <- file.path(base_path, paste0(prefix, "_", sanitized, "_count.rds"))
    meta_file <- file.path(base_path, paste0(prefix, "_", sanitized, "_meta.rds"))
    out_file <- file.path(save_path, paste0(prefix, "_", sanitized, "_cellchat.rds"))
    
    # Load count matrix and metadata from RDS files
    counts <- readRDS(count_file)
    meta <- readRDS(meta_file)
    
    # Match cell IDs between metadata and count matrix
    common <- intersect(rownames(meta), colnames(counts))
    counts <- counts[, common]; meta <- meta[common, ]
    
    # Construct and save CellChat object
    cellchat <- createCellChat(
      object = counts,
      meta = meta, #data.frame(group = meta$cluster_group, row.names = rownames(meta)),
      group.by = "cluster_group"
    )
    saveRDS(cellchat, out_file)
    message("Saved: ", out_file)
    cellchat
  })
}
