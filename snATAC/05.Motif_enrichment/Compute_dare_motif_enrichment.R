# ============================================================================
# Function: Compute Motif Enrichment for DARE
# ============================================================================
compute_dare_motif_enrichment <- function(dare, motifMat, IDX, outpath) {
  # Load required packages
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(SummarizedExperiment)
  
  # Internal enrichment function
  computeEnrichment <- function(matches = NULL, compare = NULL, background = NULL) {
    matches <- assays(matches)$matches
    
    # Compute totals
    matchCompare <- matches[compare, , drop = FALSE]
    matchBackground <- matches[background, , drop = FALSE]
    matchCompareTotal <- Matrix::colSums(matchCompare)
    matchBackgroundTotal <- Matrix::colSums(matchBackground)
    
    # Create summary data frame
    pOut <- data.frame(
      feature = colnames(matches),
      CompareFrequency = matchCompareTotal,
      nCompare = nrow(matchCompare),
      CompareProportion = matchCompareTotal / nrow(matchCompare),
      BackgroundFrequency = matchBackgroundTotal,
      nBackground = nrow(matchBackground),
      BackgroundProportion = matchBackgroundTotal / nrow(matchBackground)
    )
    
    # Enrichment
    pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProportion
    
    # Get p-values with hypergeometric test
    pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x) {
      p <- -phyper(pOut$CompareFrequency[x] - 1,
                   pOut$BackgroundFrequency[x],
                   pOut$nBackground[x] - pOut$BackgroundFrequency[x],
                   pOut$nCompare[x],
                   lower.tail = FALSE, log.p = TRUE)
      return(p / log(10))
    }) %>% unlist() %>% round(4)
    
    # Adjusted p-values
    pOut$mlog10Padj <- pmax(pOut$mlog10p - log10(ncol(pOut)), 0)
    pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]
    
    return(pOut)
  }
  
  # Step 1: Get all unique DARE peaks
  message("Step 1: Processing DARE peaks...")
  dare_all_un <- dare %>% distinct(peaks, .keep_all = TRUE)
  dare_all_idx <- rownames(IDX[IDX$peaks %in% dare_all_un$peaks, ])
  
  # Step 2: Subset motif matrix to DARE peaks
  message("Step 2: Creating DARE subset motif matrix...")
  matches <- motifMat[dare_all_idx, ]
  
  # Create output directory
  if (!dir.exists(outpath)) {
    dir.create(outpath, recursive = TRUE)
  }
  setwd(outpath)
  
  # Initialize result lists
  motif_up <- list()
  motif_down <- list()
  
  # Step 3: Process each subtype
  message("Step 3: Processing subtypes...")
  subtypes <- unique(dare$subtype)
  
  for (f in subtypes) {
    message("  Processing subtype: ", f)
    
    # Get current subtype data
    df <- dare[dare$subtype == f, ]
    
    # Process UP-regulated peaks
    df_up <- df[df$cor_trend == "Up", ]
    df_up_peaks <- df_up$peaks
    df_up_idx <- as.integer(rownames(IDX[IDX$peaks %in% df_up_peaks, ]))
    idx_up <- as.character(df_up_idx)
    
    b <- computeEnrichment(matches, idx_up, seq_len(nrow(matches)))
    b$celltype <- f
    b$trend <- "Up"
    motif_up[[f]] <- b
    write.csv(b, paste0(f, "_motif_up.csv"), row.names = FALSE)
    
    # Process DOWN-regulated peaks
    df_down <- df[df$cor_trend == "Down", ]
    df_down_peaks <- df_down$peaks
    df_down_idx <- as.integer(rownames(IDX[IDX$peaks %in% df_down_peaks, ]))
    idx_down <- as.character(df_down_idx)
    
    a <- computeEnrichment(matches, idx_down, seq_len(nrow(matches)))
    a$celltype <- f
    a$trend <- "Down"
    motif_down[[f]] <- a
    write.csv(a, paste0(f, "_motif_down.csv"), row.names = FALSE)
  }
  
  message("Processing completed! Results saved to: ", outpath)
  return(list(up = motif_up, down = motif_down))
}
