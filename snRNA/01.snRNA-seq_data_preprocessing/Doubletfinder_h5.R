library(ggplot2)
library(Seurat)
library(dplyr)
library(data.table)
library(future)
library(DoubletFinder)  # 确保已安装doubletFinder包
options(future.globals.maxSize = 50000*1024^2)
# 设置并行计算
plan(multicore, workers = 80)

# 获取当前路径下所有以"_cb_filtered_seurat.h5"结尾的h5文件
h5_files <- list.files(pattern = "*_cb_filtered_seurat\\.h5$", full.names = TRUE)

# 如果没有找到文件，尝试其他可能的模式
if (length(h5_files) == 0) {
  h5_files <- list.files(pattern = "*.h5$", full.names = TRUE)
  cat("未找到以'_cb_filtered_seurat.h5'结尾的文件，但找到以下h5文件：\n")
  print(h5_files)
}

if (length(h5_files) == 0) {
  stop("当前路径下未找到任何h5文件")
}

cat("找到", length(h5_files), "个h5文件\n")

# 创建for循环处理每个文件
for (h5_file in h5_files) {
  # 提取样本名：去掉"_cb_filtered_seurat.h5"后缀
  sample_name <- gsub("_cb_filtered_seurat\\.h5$", "", basename(h5_file))
  cat("\n===================================\n")
  cat("开始处理样本:", sample_name, "\n")
  cat("文件路径:", h5_file, "\n")
  cat("===================================\n")
  
  tryCatch({
    # 1. 读取h5文件
    cat("正在读取h5文件...\n")
    data <- Read10X_h5(h5_file)
    cat("读取完成。数据维度:", dim(data)[1], "基因 x", dim(data)[2], "细胞\n")
    
    # 2. 创建Seurat对象
    cat("创建Seurat对象...\n")
    seurat_obj <- CreateSeuratObject(
      counts = data, 
      min.cells = 3, 
      min.features = 0,
      project = sample_name
    )
    
    # 3. 添加样本信息
    seurat_obj$Sample <- sample_name
    
    # 4. 质量过滤
    cat("过滤低质量细胞 (nCount_RNA >= 500 & nFeature_RNA >= 1000)...\n")
    seurat_obj <- subset(seurat_obj, subset = nCount_RNA >= 500 & nFeature_RNA >= 1000)
    cat("过滤后细胞数:", ncol(seurat_obj), "\n")
    
    # 5. 标准化
    cat("标准化数据...\n")
    seurat_obj <- NormalizeData(
      object = seurat_obj, 
      normalization.method = "LogNormalize", 
      scale.factor = 1e4
    )
    
    # 6. 查找高变基因
    cat("查找高变基因...\n")
    seurat_obj <- FindVariableFeatures(object = seurat_obj)
    
    # 7. 缩放数据
    cat("缩放数据...\n")
    seurat_obj <- ScaleData(
      object = seurat_obj, 
      features = rownames(x = seurat_obj), 
      vars.to.regress = c("nCount_RNA")
    )
    
    # 8. PCA降维
    cat("进行PCA降维...\n")
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    
    # 9. UMAP降维
    cat("进行UMAP降维...\n")
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE, reduction = "pca")
    

    # 11. 双细胞检测
    cat("进行双细胞检测...\n")
    sweep.res.list <- paramSweep(seurat_obj, PCs = 1:10, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    
    annotations <- seurat_obj@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(0.075 * ncol(seurat_obj@assays$RNA))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    
    seurat_obj <- doubletFinder(
      seurat_obj, 
      PCs = 1:20, 
      pN = 0.25, 
      pK = mpK, 
      nExp = nExp_poi, 
      reuse.pANN = NULL, 
      sct = FALSE
    )
    
    # 12. 保存双细胞检测结果
    colnames(seurat_obj@meta.data)[length(seurat_obj@meta.data)] <- "DoubleScore"
    doublet <- table(seurat_obj$DoubleScore)  
    table_name <- paste0(sample_name, "_doublet_cell_number.txt")
    write.table(doublet, table_name, quote = FALSE, sep = "\t", row.names = FALSE)
    cat("双细胞检测结果已保存到:", table_name, "\n")
    
    # 13. 过滤双细胞
    
    seurat_obj <- subset(seurat_obj,cells=WhichCells(seurat_obj,expression=`DoubleScore`!="Doublet"))

    # 14. 保存结果
    output_file <- paste0(sample_name, "_singlet.rds")
    saveRDS(seurat_obj, output_file)
    cat("处理完成！结果已保存为:", output_file, "\n")
    
    # 15. 清理内存
    rm(data, seurat_obj)
    gc()
    
  }, error = function(e) {
    cat("处理文件时出错:", h5_file, "\n")
    cat("错误信息:", e$message, "\n")
    cat("跳过此文件，继续处理下一个...\n")
  })
}

cat("\n===================================\n")
cat("所有文件处理完成！\n")
cat("===================================\n")