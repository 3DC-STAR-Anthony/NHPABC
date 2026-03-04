#library(CellChat)
#library(Matrix)

# 处理所有年龄组的核心函数
BuildAgegroupCellChatFromCountAndMeta <- function(base_path, prefix, save_path) {
  # 年龄组和路径设置
  age_groups <- c("Old", "Exceptionally old", "Middle age", "Young")
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
  # 循环处理每个组
  lapply(age_groups, function(group) {
    # 生成文件名
    sanitized <- gsub(" ", "_", group)
    count_file <- file.path(base_path, paste0(prefix, "_", sanitized, "_count.rds"))
    meta_file <- file.path(base_path, paste0(prefix, "_", sanitized, "_meta.rds"))
    out_file <- file.path(save_path, paste0(prefix, "_", sanitized, "_cellchat.rds"))
    
    # 加载数据
    counts <- readRDS(count_file)
    meta <- readRDS(meta_file)
    
    # 匹配细胞ID
    common <- intersect(rownames(meta), colnames(counts))
    counts <- counts[, common]; meta <- meta[common, ]
    
    # 构建并保存CellChat对象
    cellchat <- createCellChat(
      object = counts,
      meta = meta, #data.frame(group = meta$cluster_group, row.names = rownames(meta)),
      group.by = "cluster_group"
    )
    saveRDS(cellchat, out_file)
    message("已保存: ", out_file)
    cellchat
  })
}
