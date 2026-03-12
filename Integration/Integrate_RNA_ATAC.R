ATACflow = function(ATAC){
    PM_ATAC_filtered <- ATAC
    
    PM_ATAC_filtered <- addClusters(
        input = PM_ATAC_filtered,
        reducedDims = "Harmony",
        method = "Seurat",
        name = "Clusters",
        resolution = 0.8,
        maxClusters = 30,
        force = TRUE
    )

    PM_ATAC_filtered <- addUMAP(
        ArchRProj = PM_ATAC_filtered, 
        reducedDims = "Harmony", 
        name = "UMAPHarmony", 
        nNeighbors = 30,       
        minDist = 0.5,         
        metric = "cosine",
        force = TRUE
    )

    PM_ATAC_filtered <- addImputeWeights(
      ArchRProj = PM_ATAC_filtered,
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
    return(PM_ATAC_filtered)
}

IntegrateFlow = function(RNA,ATAC){
    Seurat_rna <- RNA
    ArchR_atac <- ATAC
    atac_cells <- getCellNames(ArchR_atac)
    genesUse <- VariableFeatures(Seurat_rna, assay = "SCT")
    length(genesUse)
    values_to_remove <- c("TTC29", "TWSG1", "DLD", "CATSPERD", "CATSPERE")

    genesUse <- setdiff(genesUse, values_to_remove)
    useMatrix <- "GeneScoreMatrix"
    geneDF <- ArchR:::.getFeatureDF(getArrowFiles(ArchR_atac), useMatrix)
    GeneScoreMatrix <- ArchR:::.getPartialMatrix(
      ArrowFiles = getArrowFiles(ArchR_atac),
      featureDF = geneDF[geneDF$name %in% genesUse,],
      threads = 40,
      cellNames = atac_cells,
      useMatrix = useMatrix,
      verbose = FALSE
    )

    rownames(GeneScoreMatrix) <- geneDF[geneDF$name %in% genesUse, "name"]
    mat <- log(GeneScoreMatrix + 1)

    Seurat_ATAC <- CreateSeuratObject(
        counts = GeneScoreMatrix,
        assay = 'GeneScore',
        project = 'ATAC',
        min.cells = 1,
        meta.data = as.data.frame(ArchR_atac@cellColData)
    )

    Seurat_ATAC <- NormalizeData(Seurat_ATAC,verbose = FALSE)
    Seurat_ATAC <- ScaleData(Seurat_ATAC, verbose = FALSE)
    DefaultAssay(Seurat_ATAC) <- 'GeneScore'
    Seurat_ATAC$Tech <- "ATAC"
    Seurat_rna$Tech <- "RNA"

    plan("multicore", workers = 30)
    options(future.globals.maxSize = 50 * 1024^1024*1024)
    DefaultAssay(Seurat_rna) <- 'RNA'

    reference.list <- c(Seurat_rna, Seurat_ATAC)
    names(reference.list) <- c("RNA", "ATAC")
    rna_atac.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
    rna_atac_integrated <- IntegrateData(anchorset = rna_atac.anchors, dims = 1:30)
    rna_atac_integrated <- ScaleData(object = rna_atac_integrated, verbose = F)
    rna_atac_integrated <- RunPCA(object = rna_atac_integrated, verbose = F)
    rna_atac_integrated <- FindNeighbors(object = rna_atac_integrated, dims = 1:30)
    rna_atac_integrated <- FindClusters(object = rna_atac_integrated, resolution = 3)
    rna_atac_integrated <- RunHarmony(rna_atac_integrated, "Tech")
    rna_atac_integrated <- RunUMAP(rna_atac_integrated, reduction =  "harmony",reduction.name = "UMAPHarmony",dims = 1:30)
    return(rna_atac_integrated)
}
