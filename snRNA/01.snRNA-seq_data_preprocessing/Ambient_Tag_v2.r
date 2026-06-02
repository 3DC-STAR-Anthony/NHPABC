#!/usr/bin/env Rscript

# 加载所需的库
library(dplyr)
library(Seurat)
library(tidyr)
library(tidyverse)
library(DropletQC)
library(patchwork)
library(magrittr)
library(parallel)
library(GeneOverlap)
library(utils)
library(ggpubr)
source('/hwfssz1/CS_LAI/xiangyuchen/project/monkey_brain/code/AmbientRNA/FUNCTIONS.R')

# 设置全局参数
options(repr.plot.width = 12, repr.plot.height = 12)
timesExcess = 2
res = 0.5
npcs = 30
ncores = 8
ambCutoffCoef = 1.5
logfc = 0.25

# 从命令行参数获取mb参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript process_mb.R <mb>")
}
mb <- args[1]

# 拼接完整的路径
target_dir <- file.path("/hwfssz3/PS_JLU/zhangxiao6/ambient/input/new", mb)
results_dir <- file.path("/hwfssz3/PS_JLU/zhangxiao6/ambient/output/new", mb)

# 设置annotation_path
annotation_path = "/hwfssz1/CS_LAI/zhangxiao6/MBA/02.genome/ncbi_dataset/data/GCF_037993035.1/genomic.gtf"

version <- "v2"
txt_dir <- paste0(mb,"_count.csv")

sub_dirs <- c('M4_Str_1','M4_Str_2','M5_Str_1','M5_Str_2')

for (sub_dir in sub_dirs) {
  
  cat("Processing:", sub_dir, "\n")
  raw_data_dir <- file.path(target_dir,sub_dir)
  res_sub_dir <- file.path(results_dir, sub_dir)
  dir.create(res_sub_dir)
  figure_output_dir <- file.path(res_sub_dir,"Figure")
  dir.create(figure_output_dir)
  
  #STEP 1
  raw_data <- Read10X(data.dir = file.path(raw_data_dir, "FilterMatrix"), gene.column = 1)
  sampObj <- CreateSeuratObject(counts = raw_data)
  num_rows1 <- ncol(sampObj)
     
  ###STEP 2 and 3 cluster
  sampObj = subset(sampObj, subset = nCount_RNA > 0)
  finalBarcs = read.table(gzfile(file.path(raw_data_dir, "FilterMatrix/barcodes.tsv.gz")),row.names=1)
  # translate_file <- read.table(
  #                            file.path("/hwfssz3/PS_JLU/huyuebei/MonBrain/new_RNA/barcodeTranslate",paste0(sub_dir,"_barcodeTranslate_16.txt")), 
  #                            header = FALSE, 
  #                            sep = "\t", 
  #                            stringsAsFactors = FALSE)
  # colnames(translate_file) <- c("final", "filter")
  # mapping <- setNames(translate_file$final, translate_file$filter)
  # new_row_names <- mapping[rownames(finalBarcs)]
  # unique_indices <- !duplicated(new_row_names)
  # finalBarcs <- finalBarcs[unique_indices, ]
  # new_row_names <- new_row_names[unique_indices]
  # new_rownames <- unname(new_row_names[rownames(finalBarcs)])
  # new_rownames <- lapply(new_rownames, function(x) replace(x, is.na(x), 1))
  # rownames(finalBarcs) <- new_rownames
  # finalBarcs <- finalBarcs[rownames(finalBarcs) != "1", ]
                         
  sampObj$is_empty = ifelse(colnames(sampObj) %in% rownames(finalBarcs), 'Non_empty', 'Empty')
    
  if(prop.table(table(sampObj$is_empty))['Empty'] < 0.5){
		print('More than half of the cell barcodes were not filtered')
		print('We will keep all cell barcodes to find the ambient clusters')
		softFiltObj = sampObj
  } else{
		# Keep n times more than filtered
		rawMeta = sampObj@meta.data
		rawMetaFinalCount <- sum(rawMeta$is_empty == 'Non_empty')
        keepMax <- min(rawMetaFinalCount * timesExcess, nrow(rawMeta))
        keepL <- rawMeta[order(rawMeta$nCount_RNA, decreasing = TRUE)[1:keepMax], ] %>% rownames
		# Soft filter
		softFiltObj = subset(sampObj, cells = keepL)
  }
  
  
  print('Running Clustering With Excess Barcodes')
  softFiltObj = NormalizeData(softFiltObj,verbose = FALSE)
  softFiltObj = FindVariableFeatures(softFiltObj, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  softFiltObj = ScaleData(softFiltObj, features = rownames(softFiltObj),verbose = FALSE)
  softFiltObj = RunPCA(softFiltObj, verbose = FALSE)
  softFiltObj = FindNeighbors(softFiltObj, dims = 1:npcs, verbose = FALSE, reduction = 'pca')
  softFiltObj = FindClusters(softFiltObj, verbose = FALSE, resolution = res)

  filtObjMeta = softFiltObj[[]]
  clDF = filtObjMeta$seurat_clusters %>% table %>% names %>% data.frame(clusters = .)
  clDF$size = filtObjMeta[,'seurat_clusters'] %>% table
  clDF$emptySize = filtObjMeta[filtObjMeta$is_empty == 'Empty', 'seurat_clusters'] %>% table
  clDF$ratio = clDF$emptySize / clDF$size

  cutoff = prop.table(table(sampObj$is_empty))['Empty'] * ambCutoffCoef
  clDF_Ambient = clDF[clDF$ratio > cutoff | clDF$size == max(clDF$size),]

  filtObjMeta$is_ambient_cluster = 'Non_Ambient'
  filtObjMeta[filtObjMeta$seurat_clusters %in% clDF_Ambient$clusters, 'is_ambient_cluster'] = 'Ambient'
  softFiltObj@meta.data = filtObjMeta
  sampObjAmb <- softFiltObj
  sampObjAmb = RunUMAP(sampObjAmb, dims = 1:30,verbose = FALSE)
    
  ###save fig
  pdf_file <- file.path(figure_output_dir, "beforeCB_ambCluster.pdf")
  pdf(pdf_file, width = 12)
  p <- DimPlot(sampObjAmb, group.by = 'is_ambient_cluster', label = T, raster = T, label.size = 7) +theme(text=element_text(size=20, face = 'bold')) + NoLegend() + ggtitle('')
  print(p)
  dev.off()
  
  pdf_file <- file.path(figure_output_dir, "beforeCB_DotPlot.pdf")
  pdf(pdf_file, width = 12)
  p <- DotPlot(sampObjAmb, features = c('SNAP25', 'SLC17A7', 'GAD1', 'GAD2', 'PLP1', 'MBP', 'SLC1A3', 'GFAP', 'CSF1R', 'A2M', 'PDGFRA')) +theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  dev.off()
  
 
    
    
  ###STEP 5 run cellbender
  ###STEP 6 ANNOTATE MAJOR CELL TYPES AFTER CELLBENDER
  file_path <- file.path(raw_data_dir,"cb_filtered_seurat.h5")
  hmat = Read10X_h5(filename = file_path,use.names=TRUE)
  seurat_obj = CreateSeuratObject(counts = hmat)
  sample <- sub_dir %>% strsplit("_") %>% unlist() %>% .[1]
  region <- sub_dir %>% strsplit("_") %>% unlist() %>% .[2]
  library <- sub_dir %>% strsplit("_") %>% unlist() %>% .[3]
  seurat_obj <- AddMetaData(seurat_obj, metadata = data.frame(Sample = rep(sample, nrow(seurat_obj@meta.data)),
                                                              Library = rep(library, nrow(seurat_obj@meta.data)), 
                                                              Region = rep(region, nrow(seurat_obj@meta.data)),
                                                              batch = rep(version,  nrow(seurat_obj@meta.data))))
    
  seurat_obj <- NormalizeData(seurat_obj,verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj,verbose = FALSE)
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj,features = all.genes,verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj),verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15,verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.2,verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj,dims=1:15,verbose = FALSE)
  
  raw_meta <- sampObjAmb@meta.data
  new_meta <- seurat_obj@meta.data
  common_indices <- intersect(rownames(raw_meta), rownames(new_meta))
  new_meta[common_indices, "is_ambient_cluster"] <- raw_meta[common_indices, "is_ambient_cluster"]
  new_column <- new_meta$is_ambient_cluster
  seurat_obj <- AddMetaData(seurat_obj, metadata = new_column, col.name = "is_ambient_cluster")
  
  ###save fig
  pdf_file <- file.path(figure_output_dir, "afterCB_ambCluster.pdf")
  pdf(pdf_file, width = 12)
  p <- DimPlot(seurat_obj, group.by = 'is_ambient_cluster', label = T, raster = T, label.size = 7) +theme(text=element_text(size=20, face = 'bold')) + NoLegend() + ggtitle('')
  print(p)
  dev.off()
  
  pdf_file <- file.path(figure_output_dir, "afterCB_DotPlot.pdf")
  pdf(pdf_file, width = 12)
  p <- DotPlot(seurat_obj, features = c('SNAP25', 'SLC17A7', 'GAD1', 'GAD2', 'PLP1', 'MBP', 'SLC1A3', 'GFAP', 'CSF1R', 'A2M', 'PDGFRA'))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  dev.off()

  ###save rds
  saveRDS(seurat_obj, file.path(res_sub_dir, "cb_ambTag.RDS"))
  
  ###save csv
  csv_path <- file.path(results_dir,txt_dir)
  meta1 <- sampObjAmb@meta.data
  ambient_counts1 <- table(meta1$is_ambient_cluster)
  meta2 <- seurat_obj@meta.data
  num_rows2 <- nrow(meta2)
  ambient_counts2 <- table(meta2$is_ambient_cluster)
  if (file.exists(csv_path)) {
    existing_data <- read.csv(csv_path, row.names = 1)
  } else {
    existing_data <- data.frame(
    bfCB_cells = integer(),
    bfCB_Ambient = integer(),
    bfCB_Non_Ambient = integer(),
    afCB_cells = integer(),
    afCB_Ambient = integer(),
    afCB_Non_Ambient = integer(),
    stringsAsFactors = FALSE
    )
  }

  new_sample <- data.frame(
    bfCB_cells = num_rows1,
    bfCB_Ambient = ambient_counts1["Ambient"],
    bfCB_Non_Ambient = ambient_counts1["Non_Ambient"],
    afCB_cells = num_rows2,
    afCB_Ambient = ambient_counts2["Ambient"],
    afCB_Non_Ambient = ambient_counts2["Non_Ambient"],
    stringsAsFactors = FALSE
  )
  rownames(new_sample) <- sub_dir
  existing_data <- rbind(existing_data, new_sample)
  write.csv(existing_data, file = csv_path, row.names = TRUE)
    
  cat("Process End:", sub_dir, "\n")
}
