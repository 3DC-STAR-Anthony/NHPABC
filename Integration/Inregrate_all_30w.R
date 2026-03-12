#!/usr/bin/env Rscript

library(ArchR)
library(Seurat)
library(Signac)
library(future)
library(future.apply)
library(parallel)
library(dplyr)
library(harmony)
plan("multisession", workers = 1)  
options(future.globals.maxSize = 50 * 1024^1024*1024)
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 1) {
#   stop("Usage: Rscript process_mb.R <mb>")
# }
# region <- args[1]

rna_dir <- "~/RNA_downsample.rds"
atac_dir <- "~/ATAC_downsample_10w.rds"
Seurat_rna <- readRDS(rna_dir)
ArchR_atac <- readRDS(atac_dir)

Seurat_rna <- NormalizeData(Seurat_rna,verbose = FALSE)
Seurat_rna <- SCTransform(Seurat_rna, vst.flavor = "v2", verbose = TRUE, variable.features.n = 3000)
outloc1_gene <- Seurat_rna@assays$SCT@var.features[grep(invert = T,'LOC1', Seurat_rna@assays$SCT@var.features)]
Seurat_rna@assays$SCT@var.features <- unique(c(outloc1_gene))
genesUse <- VariableFeatures(Seurat_rna, assay = "SCT")
Seurat_rna <- RunPCA(Seurat_rna, features = genesUse)

atac_cells <- getCellNames(ArchR_atac)
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

plan("multicore", workers = 10)
options(future.globals.maxSize = 50 * 1024^1024*1024)
DefaultAssay(Seurat_rna) <- 'RNA'

reference.list <- c(Seurat_rna, Seurat_ATAC)
names(reference.list) <- c("RNA", "ATAC")
rna_atac.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
rna_atac_integrated <- IntegrateData(anchorset = rna_atac.anchors, dims = 1:30)
rna_atac_integrated <- ScaleData(object = rna_atac_integrated, verbose = F)
rna_atac_integrated <- RunPCA(object = rna_atac_integrated, verbose = F)
rna_atac_integrated <- FindNeighbors(object = rna_atac_integrated, dims = 1:30)
rna_atac_integrated <- FindClusters(object = rna_atac_integrated, resolution = 0.5)

rna_atac_integrated <- RunHarmony(rna_atac_integrated, "Tech")
rna_atac_integrated <- RunUMAP(rna_atac_integrated, reduction =  "harmony",reduction.name = "UMAPHarmony",dims = 1:30)

saveRDS(rna_atac_integrated, "~/Integrate_All_30w.rds")
