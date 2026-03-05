# Batch extract CellChat communication probabilities and export to CSV files
extract_cellchat_communication <- function(work_dir, pattern = ".*_customDB_processed.rds$") {
  # Set working directory to the path containing CellChat result files
  setwd(work_dir)
  
  # Get all CellChat object files matching the specified pattern
  cellchat_files <- list.files(".", pattern = pattern, full.names = TRUE)
  
  # Iterate over each CellChat file for processing
  lapply(cellchat_files, function(file_path) {
    # Extract filename (without directory path)
    file_name <- basename(file_path)
    
    # Parse tissue and age group information from filename 
    # (e.g., extract from "PFC_Old_customDB_processed.rds")
    base_name <- sub("_customDB_processed.rds$", "", file_name)
    parts <- strsplit(base_name, "_")[[1]]
    
    # Extract tissue name and age group (handle multi-word age groups)
    if (paste(tail(parts, 2), collapse = "_") %in% c("Exceptionally_old", "Middle_age")) {
      tissue <- paste(head(parts, -2), collapse = "_")
      age <- paste(tail(parts, 2), collapse = "_")
    } else {
      age_group <- tail(parts, 1)  # Single-word age groups (Old, Young)
      tissue <- paste(head(parts, -1), collapse = "_")
      age <- age_group
    }
    
    # Load CellChat object from RDS file
    message("Processing file: ", file_name)
    cellchat <- readRDS(file_path)
    
    # Extract communication probability data
    commun_df <- subsetCommunication(cellchat, thresh = 1)
    
    # Generate output filename and path
    output_file <- paste0(tissue, "_", age, "_subsetCommunication_prob_truncatedMean.csv")
    output_path <- file.path(getwd(), output_file)
    
    # Save communication data to CSV file
    write.csv(commun_df, output_path, row.names = FALSE)
    message("Saved to: ", output_path)
  })
}
