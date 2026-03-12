# Load required packages
library(Seurat)
library(ArchR)

# ============================================================================
# Function: Find Stage-Specific DAREs
# ============================================================================
# Description: This function identifies stage-specific differentially accessible
#              regions (sDAREs) between two age groups for each cell type
# Parameters:
#   - proj_path: Path to ArchR project RDS file
#   - drop_celltypes: Vector of cell types to exclude from analysis
#   - age_groups: Named list with two elements: useg (test group) and bgdg (background group)
#   - output_dir: Directory to save results
#   - pval_cutoff: p-value cutoff for marker selection. Default = 0.05
#   - threads: Number of threads for parallel processing. Default = 60
# Returns: List of data frames containing sDARE results for each cell type
# ============================================================================
find_stage_dares <- function(proj_path,
                            drop_celltypes = NULL,
                            age_groups = list(useg = "Middle", bgdg = "Old"),
                            output_dir,
                            pval_cutoff = 0.05,
                            threads = 60) {
  
  # Set working directory
  setwd(dirname(proj_path))
  
  # Load ArchR project
  message("Loading ArchR project: ", proj_path)
  proj <- readRDS(proj_path)
  
  # Get all cell types
  label3_list <- unique(proj$subtype)
  
  # Filter out excluded cell types
  if (!is.null(drop_celltypes)) {
    label3_list_f <- label3_list[!(label3_list %in% drop_celltypes)]
    message("Excluding cell types: ", paste(drop_celltypes, collapse = ", "))
  } else {
    label3_list_f <- label3_list
  }
  
  message("Processing ", length(label3_list_f), " cell types")
  message("Age groups: ", age_groups$useg, " vs ", age_groups$bgdg)
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize list for results
  all_results <- list()
  
  # Process each cell type
  for (f in label3_list_f) {
    message("\nProcessing cell type: ", f)
    
    # Step 1: Subset cells for current cell type
    message("  Step 1: Subsetting cells...")
    proj_f <- subsetCells(
      ArchRProj = proj,
      cellNames = row.names(proj@cellColData[proj@cellColData$subtype == f, ])
    )
    
    # Step 2: Find marker features
    message("  Step 2: Finding marker features...")
    markerTest <- getMarkerFeatures(
      threads = threads,
      ArchRProj = proj_f,
      useMatrix = "PeakMatrix",
      groupBy = "Age_group",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = age_groups$useg,
      bgdGroups = age_groups$bgdg
    )
    
    # Save raw marker test results
    raw_file <- paste0(output_dir, "/", f, "_", age_groups$useg, "_vs_", age_groups$bgdg, "_sDARE_raw.rds")
    message("  Step 3: Saving raw results: ", basename(raw_file))
    saveRDS(markerTest, raw_file)
    
    # Step 4: Extract significant markers
    message("  Step 4: Extracting significant markers...")
    df <- getMarkers(
      seMarker = markerTest,
      cutOff = paste("Pval <=", pval_cutoff),
      n = NULL,
      returnGR = FALSE
    )
    
    # Step 5: Format results
    df2 <- as.data.frame(df[[1]], stringsAsFactors = FALSE)
    df2$label <- f
    df2$age_stage <- paste(age_groups$useg, "_vs_", age_groups$bgdg, sep = "")
    df2$peaks <- paste(df2$seqnames, df2$start, df2$end, sep = "_")
    
    # Save formatted results
    result_file <- paste0(output_dir, "/", f, "_", age_groups$useg, "_vs_", age_groups$bgdg, "_sDARE_pvalue", pval_cutoff, ".rds")
    message("  Step 5: Saving filtered results: ", basename(result_file))
    saveRDS(df2, result_file)
    
    # Store in list
    all_results[[f]] <- df2
    
    message("  Completed: ", f, " - Found ", nrow(df2), " significant peaks")
  }
  
  message("\nAnalysis completed for ", length(all_results), " cell types")
  message("Results saved to: ", output_dir)
  
  return(all_results)
}
