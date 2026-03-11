# ============================================================================
# Function: Create Motif-Peak Matrix
# ============================================================================
# Description: This function creates a motif-peak annotation matrix from a 
#              Seurat object containing ATAC-seq data. It extracts peaks, 
#              matches motifs from JASPAR2022 database, and creates a binary
#              matrix indicating motif presence in each peak.
# Parameters:
#   - seurat_rds: Path to Seurat RDS file containing ATAC-seq data
#   - genome_annotation_rds: Path to genome annotation RDS file
#   - species: Species for motif database. Default = "Homo sapiens"
#   - collection: JASPAR collection. Default = "CORE"
#   - p_cutoff: p-value cutoff for motif matching. Default = 5e-05
#   - motif_width: Width for motif scanning. Default = 7
#   - assay_name: Name of ATAC assay in Seurat object. Default = "ATAC"
#   - slot_name: Name of data slot. Default = "data"
# Returns: A SummarizedExperiment object containing motif-peak annotation matrix
# Dependencies: Seurat, SummarizedExperiment, TFBSTools, motifmatchr, Matrix
# ============================================================================
create_motif_peak_matrix <- function(seurat_rds, 
                                     genome_annotation_rds,
                                     species = "Homo sapiens",
                                     collection = "CORE",
                                     p_cutoff = 5e-05,
                                     motif_width = 7,
                                     assay_name = "ATAC",
                                     slot_name = "data") {
  
  # Log start
  message("Creating motif-peak matrix...")
  message("  Input Seurat RDS: ", seurat_rds)
  message("  Genome annotation: ", genome_annotation_rds)
  message("  Species: ", species, ", Collection: ", collection)
  message("  Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  # Load required packages
  library(data.table)
  library(dplyr)
  library(Seurat)
  library(parallel)
  library(chromVAR)
  library(motifmatchr)
  library(SummarizedExperiment)
  library(Matrix)
  library(ggplot2)
  library(BiocParallel)
  library(ggrepel)
  library(TFBSTools)
  library(dplyr)
  library(tidyr)

  # Step 1: Load Seurat object and extract peak data
  message("Step 1: Loading Seurat object and extracting peaks...")
  seurat_obj <- readRDS(seurat_rds)
  peaks <- GetAssayData(seurat_obj, assay = assay_name, slot = slot_name)
  peak_names <- rownames(peaks)
  
  # Step 2: Parse peak names and create GRanges
  message("Step 2: Creating GRanges from peak coordinates...")
  seqnames <- sub("^([^-]+-[^-]+)-.*$", "\\1", peak_names)
  ranges <- sub("^[^-]+-[^-]+-(.*)$", "\\1", peak_names)
  seqnames <- gsub("-", "_", seqnames)
  
  rowRanges <- GRanges(
    seqnames,
    ranges,
    strand = rep("*", length(peak_names)),
    idx = 1:length(peak_names)
  )
  
  # Step 3: Get motifs from JASPAR2022
  message("Step 3: Retrieving motifs from JASPAR2022...")
  args <- list(species = species, collection = collection)
  motifs_raw <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, args)
  
  # Helper function to summarize motifs
  summarizeJASPARMotifs <- function(motifs) {
    motifNames <- lapply(seq_along(motifs), function(x) {
      namex <- make.names(motifs[[x]]@name)
      if (substr(namex, nchar(namex), nchar(namex)) == ".") {
        namex <- substr(namex, 1, nchar(namex) - 1)
      }
      namex <- paste0(namex, "_", x)
      namex
    }) %>% unlist()
    
    motifDF <- lapply(seq_along(motifs), function(x) {
      data.frame(
        row.names = motifNames[x],
        name = motifs[[x]]@name[[1]],
        ID = motifs[[x]]@ID,
        strand = motifs[[x]]@strand,
        symbol = ifelse(!is.null(motifs[[x]]@tags$symbol[1]), motifs[[x]]@tags$symbol[1], NA),
        family = ifelse(!is.null(motifs[[x]]@tags$family[1]), motifs[[x]]@tags$family[1], NA),
        alias = ifelse(!is.null(motifs[[x]]@tags$alias[1]), motifs[[x]]@tags$alias[1], NA),
        stringsAsFactors = FALSE
      )
    }) %>% Reduce("rbind", .) %>% DataFrame
    
    names(motifs) <- motifNames
    return(list(motifs = motifs, motifSummary = motifDF))
  }
  
  motif_obj <- summarizeJASPARMotifs(motifs_raw)
  motifs <- motif_obj$motifs
  motifSummary <- motif_obj$motifSummary
  message("  Retrieved ", length(motifs), " motifs")
  
  # Step 4: Load genome annotation
  message("Step 4: Loading genome annotation...")
  genomeAnnotation <- readRDS(genome_annotation_rds)
  
  # Step 5: Match motifs to peak regions
  message("Step 5: Matching motifs to peak regions...")
  motifPositions <- motifmatchr::matchMotifs(
    pwms = motifs,
    subject = rowRanges,
    genome = genomeAnnotation$genome,
    out = "positions",
    p.cutoff = p_cutoff,
    w = motif_width
  )
  
  # Step 6: Create motif-peak matrix
  message("Step 6: Creating motif-peak sparse matrix...")
  allPositions <- unlist(motifPositions)
  overlapMotifs <- findOverlaps(rowRanges, allPositions, ignore.strand = TRUE)
  
  motifMat <- Matrix::sparseMatrix(
    i = queryHits(overlapMotifs),
    j = match(names(allPositions), names(motifPositions))[subjectHits(overlapMotifs)],
    x = rep(TRUE, length(overlapMotifs)),
    dims = c(length(rowRanges), length(motifPositions))
  )
  colnames(motifMat) <- names(motifPositions)
  motifMat <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(matches = motifMat), rowRanges = rowRanges)
  rownames(motifMat) <- 1:length(peak_names)
  saveRDS(motif_matrix, paste0(output_path,"NHPABC_allpeak_motif_annotation_matrix.rds"))

  return(motif_matrix)
}
