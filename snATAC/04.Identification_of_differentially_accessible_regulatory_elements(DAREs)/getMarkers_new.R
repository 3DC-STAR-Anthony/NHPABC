getMarkers_new <- function(
  seMarker = NULL,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  n = NULL,
  returnGR = FALSE
  ){

  # 简单的输入验证替代 .validInput
  if (is.null(seMarker)) {
    stop("seMarker cannot be NULL")
  }
  if (!inherits(seMarker, "SummarizedExperiment")) {
    stop("seMarker must be a SummarizedExperiment object")
  }
  if (!is.character(cutOff)) {
    stop("cutOff must be a character string")
  }
  if (!is.null(n) && !is.numeric(n)) {
    stop("n must be NULL or a numeric value")
  }
  if (!is.logical(returnGR)) {
    stop("returnGR must be a logical value (TRUE or FALSE)")
  }

  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }

  if(returnGR){

    if(S4Vectors::metadata(seMarker)$Params$useMatrix != "PeakMatrix"){
      stop("Only markers can be returned as GRanges when PeakMatrix!")
    }

    rr <- GRanges(rowData(seMarker)$seqnames, IRanges(rowData(seMarker)$start, rowData(seMarker)$end))

    grL <- lapply(seq_len(ncol(passMat)), function(x){
      idx <- which(passMat[, x])
      rrx <- rr[idx]
      rrx$Log2FC <- SummarizedExperiment::assays(seMarker[idx, ])[["Log2FC"]][, x]
      rrx$FDR <- SummarizedExperiment::assays(seMarker[idx, ])[["FDR"]][, x]
      # 添加Pval提取
      if("Pval" %in% assayNames){
        rrx$Pval <- SummarizedExperiment::assays(seMarker[idx, ])[["Pval"]][, x]
      }
      if("MeanDiff" %in% assayNames){
        rrx$MeanDiff <- SummarizedExperiment::assays(seMarker[idx, ])[["MeanDiff"]][, x]
      }
      rrx <- rrx[order(rrx$FDR),,drop=FALSE]
      if(!is.null(n)){
        if(n < nrow(rrx)){
          rrx <- rrx[seq_len(n), , drop = FALSE]
        }
      }
      rrx
    }) %>% SimpleList

    names(grL) <- colnames(seMarker)

    grL <- grL[gtools::mixedsort(names(grL))]

    return(grL)

  }else{

    markerList <- lapply(seq_len(ncol(passMat)), function(x){
      idx <- which(passMat[, x])
      rrx <- SummarizedExperiment::rowData(seMarker[idx,])
      rrx$Log2FC <- SummarizedExperiment::assays(seMarker[idx, ])[["Log2FC"]][, x]
      rrx$FDR <- SummarizedExperiment::assays(seMarker[idx, ])[["FDR"]][, x]
      # 添加Pval提取
      if("Pval" %in% assayNames){
        rrx$Pval <- SummarizedExperiment::assays(seMarker[idx, ])[["Pval"]][, x]
      }
      if("MeanDiff" %in% assayNames){
        rrx$MeanDiff <- SummarizedExperiment::assays(seMarker[idx, ])[["MeanDiff"]][, x]
      }
      rrx <- rrx[order(rrx$FDR),,drop=FALSE]
      if(!is.null(n)){
        if(n < nrow(rrx)){
          rrx <- rrx[seq_len(n), , drop = FALSE]
        }
      }
      rrx
    }) %>% SimpleList

    names(markerList) <- colnames(seMarker)

    return(markerList)

  }

}
