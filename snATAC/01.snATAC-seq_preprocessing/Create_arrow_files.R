# Load required packages
library(RColorBrewer)
library(ArchR)
library(dplyr)
library(ggplot2)
library(BSgenome.Mfascicularis.NCBI.T2TMFA8v1)
library(GenomicFeatures)
library(OrganismDbi)

# ============================================================================
# Function: Create ArchR Arrow Files
# ============================================================================
# Description: This function creates ArchR Arrow files from input fragment files
# Parameters:
#   - input_dir: Directory containing input fragment files (.gz format)
#   - genome_annotation_path: Path to genome annotation RDS file
#   - gene_annotation_path: Path to gene annotation RDS file
#   - minTSS: Minimum TSS score for filtering. Default = 4
#   - minFrags: Minimum fragments per cell. Default = 3000
#   - threads: Number of threads for parallel processing. Default = 40
#   - pattern: File pattern to match. Default = "\\.fragments\\.tsv\\.gz$"
# Returns: List of Arrow file paths
# ============================================================================
create_arrow_files <- function(input_dir, 
                               genome_annotation_path, 
                               gene_annotation_path,
                               minTSS = 4,
                               minFrags = 3000,
                               threads = 40,
                               pattern = "\\.fragments\\.tsv\\.gz$") {
  
  # Set seed and working directory
  set.seed(1)
  setwd(input_dir)
  
  # Set ArchR parameters
  addArchRThreads(threads = threads)
  addArchRChrPrefix(chrPrefix = FALSE)
  addArchRLocking(locking = FALSE)
  
  # Set plot options
  options(repr.plot.width = 12, repr.plot.height = 12)
  
  # Get input files
  args <- dir(pattern = pattern)
  name <- gsub("\\.fragments\\.tsv\\.gz$", "", args)
  inputFiles <- args
  sampleNames <- name
  
  message("Found ", length(inputFiles), " input files")
  message("Sample names: ", paste(sampleNames, collapse = ", "))
  
  # Load genome and gene annotations
  message("Loading genome and gene annotations...")
  genomeAnnotation <- readRDS(genome_annotation_path)
  geneAnnotation <- readRDS(gene_annotation_path)
  
  # Create Arrow files
  message("Creating Arrow files...")
  ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = sampleNames,
    minTSS = minTSS,
    minFrags = minFrags,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
    subThreading = TRUE,
    threads = threads
  )
  
  message("Arrow files created successfully!")
  message("Arrow file paths: ", paste(ArrowFiles, collapse = ", "))
  
  return(ArrowFiles)
}
