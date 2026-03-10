#' Calculate and filter cCRE
#'
#' This function takes a Seurat object (peak matrix) for a given cell type,
#' computes mean CPM values per monkey individual, applies filtering criteria,
#' and returns the list of enriched peaks.
#'
#' @param rds_path Path to the input Seurat RDS file (peak matrix for one cell type)
#' @param output_dir Directory to save output cCRE RDS files
#' @param min_samples_cpm2 Minimum number of samples with mean CPM > 4 (default: 4)
#' @param min_samples_cpm0 Minimum number of samples with mean CPM > 0 (default: 12)
#'
#' @return A data frame of filtered cCRE values, also saved as an RDS file
#'
#' @import Seurat
#' @export
calculate_cCRE <- function(rds_path,
                          output_dir = "./03.Cell_type_enrich_peak/",
                          min_samples_cpm2 = 4,
                          min_samples_cpm0 = 12) {

  # Check if input file exists
  if (!file.exists(rds_path)) {
    stop("Input RDS file does not exist: ", rds_path)
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Extract cell type name from file path (adjust regex if needed)
  cell_type <- sub("_peakmatrix_seurat\\.rds$", "", basename(rds_path))

  # Read Seurat object
  rds <- readRDS(rds_path)
  data_matrix <- rds@assays$peaks@data
  meta_data <- rds@meta.data

  # Initialize data frame with peaks as rownames
  df <- data.frame(peaks = rownames(data_matrix))
  rownames(df) <- df$peaks
  df <- df[, -1, drop = FALSE]

  # Calculate mean CPM for each monkey individual
  for (sample in unique(meta_data$individual)) {
    # Get cells belonging to this individual
    sample_cells <- rownames(meta_data[meta_data$individual == sample, ])
    # Subset matrix to these cells
    sample_matrix <- data_matrix[, sample_cells, drop = FALSE]
    # Compute row means (mean CPM per peak for this individual)
    sample_mean_cpm <- rowMeans(sample_matrix)
    # Add to data frame
    df[[sample]] <- sample_mean_cpm
  }

  # Save intermediate mean CPM matrix
  saveRDS(df, file.path(output_dir, paste0(cell_type, "_cCRE_row.rds")))

  # Apply filtering criteria
  # Condition 1: mean CPM > 2 in at least min_samples_cpm2 samples
  peak_f1 <- df[rowSums(df > 2) >= min_samples_cpm2, , drop = FALSE]
  # Condition 2: mean CPM > 0 in at least min_samples_cpm0 samples
  peak_f2 <- df[rowSums(df > 0) >= min_samples_cpm0, , drop = FALSE]

  # Find intersect of peaks passing both filters
  filtered_peaks <- intersect(rownames(peak_f1), rownames(peak_f2))

  # Extract final CEP data frame
  cep <- df[filtered_peaks, , drop = FALSE]

  # Save final CEP result
  saveRDS(cep, file.path(output_dir, paste0(cell_type, "_cCRE.rds")))

  # Print summary
  message("Cell type: ", cell_type)
  message("Total peaks after filtering: ", length(filtered_peaks))

  # Return the CEP data frame
  return(cep)
}
