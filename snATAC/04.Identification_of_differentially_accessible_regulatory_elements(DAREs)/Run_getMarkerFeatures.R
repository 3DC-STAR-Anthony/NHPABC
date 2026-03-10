#!/usr/bin/env Rscript

library(ArchR)
library(parallel)
library(ggplot2)

# Set parallel processing
addArchRThreads(threads = 10)

# Load ArchR project
cat("Loading ArchR project...\n")
proj <- readRDS("~/ArchR_project.rds")

# Modify cell metadata: create LongDARE groups
# Mark "Exceptionally old" age group as "old", others as "rest"
proj@cellColData$LongDARE_group <- ifelse(proj@cellColData$Age_group == "Exceptionally old", "old", "rest")

# Get all cell types
df_celltype <- data.frame(celltype = sort(unique(proj$subtype)))
cell_types <- df_celltype[, "celltype"]

cat("Cell types to process:", paste(cell_types, collapse = ", "), "\n")

# Process each cell type
for (f in cell_types) {
    cat("\nProcessing cell type:", f, "------", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    
    # Step 1: Filter cells of current cell type
    idx_celltype <- which(proj$subtype == f)
    proj_subset <- proj[idx_celltype, ]
    
    # Step 2: Identify differential peaks using getMarkerFeatures
    cat("  Running getMarkerFeatures...\n")
    markers <- getMarkerFeatures(
        ArchRProj = proj_subset,           # ArchR project subset
        useMatrix = "PeakMatrix",          # Use peak matrix
        groupBy = "LongDARE_group",        # Group by LongDARE
        useGroups = "old",                 # Test group: "old" (long-lived)
        bgdGroups = "rest",                # Background group: "rest" (non-long-lived)
        testMethod = "wilcoxon",           # Use Wilcoxon rank-sum test
        normBy = "ReadsInTSS",             # Normalize by TSS reads
        bias = c("TSSEnrichment", "log10(nFrags)"),  # Correct for TSS enrichment and fragment count bias
        useSeqnames = NULL,                # Use all chromosomes
        verbose = TRUE,                    # Show detailed log
        logFile = createLogFile(paste0("getMarkerFeatures_", f))  # Log file
    )
    
    # Step 3: Save results
    output_file <- paste0("~/NHPABC/Figure6/07.LongDAREs/01.getMarkerFeatures_result/", f, ".LongDAREs_raw.rds")
    saveRDS(markers, output_file)
    cat("  Saved to:", output_file, "\n")
}

cat("\nBatch 1 completed!\n")
cat("End time:", as.character(Sys.time()), "\n")
