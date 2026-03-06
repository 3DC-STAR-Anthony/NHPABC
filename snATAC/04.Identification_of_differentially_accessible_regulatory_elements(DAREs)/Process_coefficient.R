# ============================================================================
# Function: Process Coefficient Analysis from RDS File
# ============================================================================
# Description: This function reads CPM matrix from RDS file, performs coefficient analysis.
#              It creates Seurat object, applies age grouping, converts to SingleCellExperiment,
#              runs MAST analysis, and saves results.
# Parameters:
#   - input_rds: Path to the input RDS file containing CPM matrix
#   - output_dir: Directory to save the result RDS file
# Returns: A data.frame containing coefficient analysis results
# Dependencies: Seurat, SingleCellExperiment, MAST, data.table
# ============================================================================
process_coefficient <- function(input_rds, output_dir) {
  # Extract subtype name from filename
  subtype_name <- tools::file_path_sans_ext(basename(input_rds))
  
  # Log start of processing
  message("Processing coefficient analysis: ", subtype_name, 
          " ------ ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  # Read CPM matrix from RDS
  message("  Reading CPM matrix from: ", input_rds)
  cpm_matrix <- readRDS(input_rds)
  message("  Matrix dimensions: ", nrow(cpm_matrix), " features x ", ncol(cpm_matrix), " samples")
  
  # Define sample information
  monkey <- c('Y1','Y2','Y3','Y4','Y5','Y6','Y7',
              'M1','M2','M3','M4','M5',
              'O1','O2','O3','O4','O5','O6',
              'V1','V2','V3','V4','V5','V6')
  aging <- c(5,6,5,5,5,5,5,
             10,11,12,10,10,
             22,23,22,22,22,23,
             29,29,28,31,29,31)
  monkey_age <- data.frame(monkey = monkey, age = aging)
  rownames(monkey_age) <- monkey_age$monkey
  
  # Prepare metadata
  message("  Step 1: Preparing metadata...")
  meta_s <- monkey_age[colnames(cpm_matrix), ]
  
  # Create Seurat object
  message("  Step 2: Creating Seurat object...")
  filter <- CreateSeuratObject(
    counts = as(as.matrix(cpm_matrix), "sparseMatrix"),
    assay = "peaks",
    meta.data = meta_s
  )
  
  # Group age values
  message("  Step 3: Grouping age values...")
  filter$age[filter$age == '6'] = '5'
  filter$age[filter$age == '10'] = '11'
  filter$age[filter$age == '12'] = '11'
  filter$age[filter$age == '23'] = '22'
  filter$age[filter$age == '28'] = '31'
  filter$age[filter$age == '31'] = '31'
  
  # Convert to SingleCellExperiment
  message("  Step 4: Converting to SingleCellExperiment...")
  rds.sce <- as.SingleCellExperiment(filter)
  sca <- SceToSingleCellAssay(rds.sce, class = "SingleCellAssay", check_sanity = FALSE)
  
  # Filter zero rows
  message("  Step 5: Filtering zero rows...")
  sca_filt = sca[rowSums(assay(sca)) != 0, ]
  
  # Clean up memory
  rm(rds.sce, sca, filter)
  gc()
  
  # Prepare for MAST analysis
  message("  Step 6: Preparing for MAST analysis...")
  sca_filt$age <- as.numeric(sca_filt$age)
  
  # Run MAST zero-inflated model
  message("  Step 7: Running MAST zero-inflated model...")
  zlmCond <- zlm(~age, sca = sca_filt)
  summaryCond <- summary(zlmCond, doLRT = "age")
  summaryDt <- summaryCond$datatable
  
  # Extract H component results
  message("  Step 8: Extracting H component results...")
  dt1 = summaryDt[contrast == "age" & component == "H", .(primerid, `Pr(>Chisq)`)]
  
  # Extract logFC component results
  message("  Step 9: Extracting logFC component results...")
  dt2 = summaryDt[contrast == "age" & component == "logFC", .(primerid, coef, z)]
  
  # Clean up memory
  rm(sca_filt, zlmCond, summaryCond, summaryDt)
  gc()
  
  # Merge results
  message("  Step 10: Merging results...")
  de_res = merge(dt1, dt2, by = "primerid")
  colnames(de_res) <- c("gene", "age.H_p", "age.logFC", 'age.logFC_z')
  
  # Calculate FDR
  message("  Step 11: Calculating FDR...")
  de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "fdr")
  
  # Add subtype label
  de_res$label <- subtype_name
  de_res <- data.frame(de_res)
  
  # Save results
  output_path <- file.path(output_dir, paste0(subtype_name, "_coefficient.rds"))
  message("  Step 12: Saving results to: ", output_path)
  saveRDS(de_res, output_path)
  
  # Log completion
  message("Completed processing: ", subtype_name, 
          " ------ ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  return(de_res)
}
