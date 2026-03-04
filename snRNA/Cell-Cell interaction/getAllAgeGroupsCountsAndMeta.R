getAllAgeGroupsCountsAndMeta <- function(obj, cluster, ident = "Age_group", save_path = NULL, prefix = NULL) {
  # 定义需要处理的所有Age_group分组
  allowed_ages <- c("Old", "Exceptionally old", "Middle age", "Young")
  
  # 纠错检查：必要参数检查
  if (missing(cluster)) stop("请提供'cluster'参数（metadata中细胞分组列名）")
  if (missing(save_path)) stop("请提供保存路径（save_path参数）")
  if (missing(prefix)) stop("请提供文件名前缀（prefix参数）")
  
  # 确保保存路径存在
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    message("已创建保存目录: ", save_path)
  }
  
  # 纠错检查：Seurat对象metadata存在性
  if (!inherits(obj, "Seurat")) stop("输入必须是Seurat对象")
  meta <- obj@meta.data
  if (!ident %in% colnames(meta)) {
    stop(paste("metadata中未找到", ident, "列，请检查列名拼写"))
  }
  
  # 纠错检查：cluster列是否存在
  if (!cluster %in% colnames(meta)) {
    stop(paste("metadata中未找到", cluster, "列，请检查列名拼写"))
  }
  
  # 存储所有结果的列表
  all_results <- list()
  
  # 循环处理每个年龄分组
  for (group in allowed_ages) {
    message("\n处理分组: ", group)
    
    # 核心逻辑：根据group筛选细胞
    obj_sub <- subset(obj, subset = !!sym(ident) == group)
    
    # 提取count矩阵（使用Seurat V5的layer参数）
    count_matrix <- GetAssayData(obj_sub, assay = "RNA", layer = "data") # normalized data matrix
    
    # 提取元数据并保留cluster信息
    meta_data <- obj_sub@meta.data[, c(cluster, ident), drop = FALSE]
    colnames(meta_data)[colnames(meta_data) == cluster] <- "cluster_group"
    
    # 生成标准化文件名（包含前缀）
    sanitized_group <- gsub(" ", "_", group)  # 替换空格为下划线
    count_filename <- paste0(prefix, "_", sanitized_group, "_count.rds")
    meta_filename <- paste0(prefix, "_", sanitized_group, "_meta.rds")
    
    count_path <- file.path(save_path, count_filename)
    meta_path <- file.path(save_path, meta_filename)
    
    # 保存count矩阵（以RDS格式保存稀疏矩阵）
    saveRDS(count_matrix, file = count_path)
    
    # 保存元数据（以RDS格式保存，保留行名和数据类型）
    saveRDS(meta_data, file = meta_path)
    
    # 打印保存信息
    message("count矩阵已保存至: ", count_path)
    message("元数据已保存至: ", meta_path)
    
    # 存储结果到列表
    all_results[[group]] <- list(
      count_matrix = count_matrix,
      meta_data = meta_data#,
      #count_path = count_path,
      #meta_path = meta_path
    )
  }
  
  message("\n所有分组处理完成")
  return(all_results)
}

# 使用示例：
# 假设已经有一个名为seurat_obj的Seurat V5对象
# getAllAgeGroupsCountsAndMeta(
#   obj = seurat_obj,
#   cluster = "cell_type",  # 替换为实际的细胞分组列名
#   ident = "Age_group",    # 年龄分组列名，默认即可
#   save_path = "./output", # 保存目录
#   prefix = "liver"        # 文件名前缀
# )
# 会生成以下文件:
# ./output/liver_Old_count.rds
# ./output/liver_Old_meta.rds
# ./output/liver_Exceptionally_old_count.rds
# ./output/liver_Exceptionally_old_meta.rds
# ./output/liver_Middle_age_count.rds
# ./output/liver_Middle_age_meta.rds
# ./output/liver_Young_count.rds
# ./output/liver_Young_meta.rds
    
