getAllAgeGroupsCountsAndMeta <- function(obj, cluster, ident = "Age_group", save_path = NULL, prefix = NULL) {
  # Define all Age_group categories to be processed
  allowed_ages <- c("Old", "Exceptionally old", "Middle age", "Young")
  
  # Error checking: Mandatory parameter validation
  if (missing(cluster)) stop("Please provide the 'cluster' parameter (column name of cell grouping in metadata)")
  if (missing(save_path)) stop("Please provide the save path (save_path parameter)")
  if (missing(prefix)) stop("Please provide the filename prefix (prefix parameter)")
  
  # Ensure save directory exists
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    message("Created save directory: ", save_path)
  }
  
  # Error checking: Validate input is a Seurat object and metadata columns exist
  if (!inherits(obj, "Seurat")) stop("Input must be a Seurat object")
  meta <- obj@meta.data
  if (!ident %in% colnames(meta)) {
    stop(paste("Column", ident, "not found in metadata. Please check column name spelling"))
  }
  
  # Error checking: Validate cluster column exists in metadata
  if (!cluster %in% colnames(meta)) {
    stop(paste("Column", cluster, "not found in metadata. Please check column name spelling"))
  }
  
  # List to store all processing results
  all_results <- list()
  
  # Iterate through each age group for processing
  for (group in allowed_ages) {
    message("\nProcessing group: ", group)
    
    # Core logic: Subset Seurat object by age group
    obj_sub <- subset(obj, subset = !!sym(ident) == group)
    
    # Extract count matrix (using Seurat V5 layer parameter)
    count_matrix <- GetAssayData(obj_sub, assay = "RNA", layer = "data") # normalized data matrix
    
    # Extract metadata and retain cluster information
    meta_data <- obj_sub@meta.data[, c(cluster, ident), drop = FALSE]
    colnames(meta_data)[colnames(meta_data) == cluster] <- "cluster_group"
    
    # Generate sanitized filenames (with prefix)
    sanitized_group <- gsub(" ", "_", group)  # Replace spaces with underscores
    count_filename <- paste0(prefix, "_", sanitized_group, "_count.rds")
    meta_filename <- paste0(prefix, "_", sanitized_group, "_meta.rds")
    
    count_path <- file.path(save_path, count_filename)
    meta_path <- file.path(save_path, meta_filename)
    
    # Save count matrix (save sparse matrix in RDS format)
    saveRDS(count_matrix, file = count_path)
    
    # Save metadata (save in RDS format, retain row names and data types)
    saveRDS(meta_data, file = meta_path)
    
    # Print save information
    message("Count matrix saved to: ", count_path)
    message("Metadata saved to: ", meta_path)
    
    # Store results in list
    all_results[[group]] <- list(
      count_matrix = count_matrix,
      meta_data = meta_data#,
      #count_path = count_path,
      #meta_path = meta_path
    )
  }
  
  message("\nAll groups processed successfully")
  return(all_results)
}

# Usage example:
# Assume you have a Seurat V5 object named seurat_obj
# getAllAgeGroupsCountsAndMeta(
#   obj = seurat_obj,
#   cluster = "cell_type",  # Replace with actual cell grouping column name
#   ident = "Age_group",    # Age grouping column name (default is acceptable)
#   save_path = "./output", # Save directory
#   prefix = "liver"        # Filename prefix
# )
# This will generate the following files:
# ./output/liver_Old_count.rds
# ./output/liver_Old_meta.rds
# ./output/liver_Exceptionally_old_count.rds
# ./output/liver_Exceptionally_old_meta.rds
# ./output/liver_Middle_age_count.rds
# ./output/liver_Middle_age_meta.rds
# ./output/liver_Young_count.rds
# ./output/liver_Young_meta.rds
