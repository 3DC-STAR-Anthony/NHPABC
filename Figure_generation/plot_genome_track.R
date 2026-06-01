#Full Script Description
#This script does NOT directly plot genome tracks. It only contains code to generate bigWig (.bw) files from an ArchR project for genome browser visualization.
#Full visualization workflow:
#Step1: Run this script to generate bw files by cell type × age group from your ArchR project.
#Step2: Open IGV genome browser and load the matching reference genome (hg38 / mm10).
#Step3: Load all bw files into IGV and navigate to your target gene / peak region.
#Step4: Adjust track colors and signal scales in IGV, then export high-quality images for figures.
#Grouping logic:
#The ArchR project contains two metadata columns: celltype (cell type) and age_group (age group).
#These columns are combined to create cell-type-specific and age-specific pseudo-bulk groups for bw generation.

# Load required package
library(ArchR)

# Basic settings
addArchRGenome("hg38")
addArchRThreads(threads = c(8))

# Load ArchR project
proj <- readRDS("./ArchR_Project.rds")

# Create group: cell type + age group
proj$cell_age <- paste0(proj$celltype, "_", proj$age_group)

# Filter low-cell groups
min_cells <- c(50)
valid_groups <- names(table(proj$cell_age))[table(proj$cell_age) >= min_cells]
proj <- proj[proj$cell_age %in% valid_groups, ]

# Add group coverages (required for bigWig)
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "cell_age",
  minCells = c(50),
  maxCells = c(2000)
)

# Generate normalized bigWig files
getGroupBW(
  ArchRProj = proj,
  groupBy = "cell_age",
  normMethod = "ReadsInTSS",
  tileSize = c(100),
  ceiling = c(5)
)

# Export bw file list for IGV
bw_dir <- paste0(proj@projectMetadata$outputDirectory, "/GroupBigWigs/cell_age/")
bw_files <- list.files(bw_dir, ".bw$", full.names = TRUE)
writeLines(bw_files, "bw_file_list.txt")

# Message
message("BigWig files generated! Load into IGV for genome track visualization.")
