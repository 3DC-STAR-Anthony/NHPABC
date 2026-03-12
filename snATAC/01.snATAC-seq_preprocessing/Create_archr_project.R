# ============================================================================
# Function: Create ArchR Project
# ============================================================================
# Description: This function creates an ArchR project from fragment files,
#              performs quality control, clustering, and dimensionality reduction
# Parameters:
#   - input_dir: Directory containing input fragment files (.gz format)
#   - genome_annotation_path: Path to genome annotation RDS file
#   - gene_annotation_path: Path to gene annotation RDS file
#   - output_dir: Output directory for saving results
#   - minTSS: Minimum TSS score for filtering. Default = 4
#   - minFrags: Minimum fragments per cell. Default = 3000
#   - threads: Number of threads for parallel processing. Default = 20
# Returns: ArchR project object
# ============================================================================
create_archr_project <- function(input_dir, 
                                 genome_annotation_path, 
                                 gene_annotation_path,
                                 output_dir,
                                 minTSS = 4,
                                 minFrags = 3000,
                                 threads = 20) {
  
  # Load required packages
  library(RColorBrewer)
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(BSgenome.Mfascicularis.NCBI.T2TMFA8v1)
  library(GenomicFeatures)
  library(OrganismDbi)
  
  # Set seed and working directory
  set.seed(1)
  setwd(input_dir)
  
  # Set ArchR parameters
  addArchRThreads(threads = threads)
  addArchRChrPrefix(chrPrefix = FALSE)
  addArchRLocking(locking = FALSE)
  
  # Step 1: Find input files
  message("Step 1: Finding input fragment files...")
  args <- dir(pattern = "\\.gz$")
  name <- gsub("\\.fragments\\.tsv\\.gz$", "", args)
  inputFiles <- args
  sampleNames <- name
  
  message("  Input files: ", length(inputFiles))
  message("  Sample names: ", paste(sampleNames, collapse = ", "))
  
  # Load genome and gene annotations
  message("Step 2: Loading genome and gene annotations...")
  genomeAnnotation <- readRDS(genome_annotation_path)
  geneAnnotation <- readRDS(gene_annotation_path)
  
  # Step 3: Create Arrow files
  message("Step 3: Creating Arrow files...")
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
  
  # Step 4: Add doublet scores
  message("Step 4: Adding doublet scores...")
  doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10,
    knnMethod = "UMAP",
    LSIMethod = 1
  )
  
  # Step 5: Create ArchR project
  message("Step 5: Creating ArchR project...")
  proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = output_dir,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
    threads = threads,
    copyArrows = TRUE
  )
  
  # Step 6: Filter doublets
  message("Step 6: Filtering doublets...")
  proj <- filterDoublets(proj, filterRatio = 2)
  
  # Step 7: Perform iterative LSI
  message("Step 7: Performing iterative LSI...")
  proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 5,
    clusterParams = list(
      resolution = c(0.2),
      sampleCells = 10000,
      n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30
  )
  
  # Step 8: Apply Harmony batch correction
  message("Step 8: Applying Harmony batch correction...")
  proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    theta = 10,
    lambda = 10,
    nclust = 100,
    max_iter = 20,
    force = TRUE
  )
  
  # Step 9: Add clusters
  message("Step 9: Adding clusters...")
  proj <- addClusters(
    input = proj,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    maxClusters = 30,
    force = TRUE
  )
  
  # Step 10: Add UMAP
  message("Step 10: Adding UMAP visualization...")
  proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
  )
  
  # Step 11: Add imputation weights
  message("Step 11: Adding imputation weights...")
  proj <- addImputeWeights(
    ArchRProj = proj,
    reducedDims = "Harmony",
    dimsToUse = NULL,
    scaleDims = NULL,
    corCutOff = 0.75,
    td = 3,
    ka = 4,
    sampleCells = 5000,
    nRep = 2,
    k = 15,
    epsilon = 1,
    useHdf5 = TRUE,
    randomSuffix = FALSE,
    threads = getArchRThreads(),
    seed = 1,
    verbose = TRUE,
    logFile = createLogFile("addImputeWeights")
  )
  
  # Step 12: Save project
  message("Step 12: Saving project to RDS file...")
  saveRDS(proj, paste0(output_dir, "/ArchR_project.rds"))
  message("Project saved to: ", paste0(output_dir, "/ArchR_project.rds"))
  
  message("ArchR project creation completed successfully!")
  return(proj)
}
