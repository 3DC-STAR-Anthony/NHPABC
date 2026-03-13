library(Seurat)
library(tidyverse)
library(dplyr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggsci)
library(parallel)
library(RUVSeq)
library(Seurat)
library(DESeq2)
library(S4Vectors)
library(stats4)
library(BiocGenerics)
library(TRADEtools)

rds <- readRDS('/hwfssz3/PS_JLU/zhangxiao6/scRNA/Each_single_region/')


protein_cd <- read.csv('/hwfssz3/PS_JLU/zhangxiao6/scRNA/protein_coding_gene_list.csv')

anno_gene <- unique(protein_cd$Gene)[grep('LOC1',unique(protein_cd$Gene),invert = T)]
filter <- subset(rds, features = anno_gene)
deg_list <- read.csv('/hwfssz3/PS_JLU/zhangxiao6/Figure3_DEG/Revise_DEG_LIST/Revised_DEG_all_filter.csv')



# 定义要循环的参数
#deg_list <- read.csv('/hwfssz3/PS_JLU/zhangxiao6/Figure3_DEG/Revise_DEG_LIST/DEG_ALL_filter_new_opc.csv')
outdir <- '/hwfssz3/PS_JLU/zhangxiao6/Figure3_DEG/Trade2/Trade_DEseq2/'
outdir1 <- '/hwfssz3/PS_JLU/zhangxiao6/Figure3_DEG/Trade2/Trade_result_list/'
outdir2 <- '/hwfssz3/PS_JLU/zhangxiao6/Figure3_DEG/Trade2/Trade_result_table/'
#deg_list <- read.csv('/hwfssz3/PS_JLU/zhangxiao6/Figure3_DEG/Revise_DEG_LIST/DEG_ALL_filter_new_opc.csv')
Region_list <- unique(rds$Region)
subtypes <- unique(deg_list$subtype[deg_list$Region == Region_list])
age_combinations <- list(
  c('Young', 'Middle age'),
  c('Middle age', 'Old'),
  c('Old', 'Exceptionally old')
)
# 函数：清理字符串中的非法字符
clean_filename <- function(name) {
  # 替换特殊字符
  clean_name <- gsub('[+/\\\\:*?"<>| ]', '_', name)
  # 去除首尾的下划线
  clean_name <- gsub('^_+|_+$', '', clean_name)
  # 将多个连续下划线替换为单个
  clean_name <- gsub('_+', '_', clean_name)
  return(clean_name)
}

# 循环处理每个subtype和年龄组合
for(f in subtypes) {
  for(age_pair in age_combinations) {
    Age_group1 <- age_pair[1]
    Age_group2 <- age_pair[2]
    
    # 您的原始代码（稍作修改）
    filter1 <- subset(filter, subtype == f)
    filter2 <- subset(filter1, Age_group == Age_group1 | Age_group == Age_group2)
    
    # 检查是否有足够的数据
    if(nrow(filter2@meta.data) == 0) {
      message(paste("跳过:", f, Age_group1, "vs", Age_group2, "- 无数据"))
      next
    }
    
    aggr <- AggregateExpression(filter2, 
                                assays = 'RNA',
                                group.by = c("Age_group","Sample"),
                                return.seurat = FALSE)
    
    # 检查是否有足够的样本进行DESeq2分析
    if(length(unique(filter2$Age_group)) < 2) {
      message(paste("跳过:", f, Age_group1, "vs", Age_group2, "- 年龄组不足"))
      next
    }
    
    av <- as.data.frame(aggr[[1]])
    countData <- av
    
    # 检查是否有足够的列
    if(ncol(countData) < 2) {
      message(paste("跳过:", f, Age_group1, "vs", Age_group2, "- 样本数不足"))
      next
    }
    
    colData <- str_split_fixed(colnames(countData), '_', n = 2) %>% 
      as.data.frame() %>% 
      setNames(c('Age_group', 'Sample'))
    colData$id <- colnames(countData)
    rownames(colData) <- colData$id
    
    # DESeq2分析
    tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = round(countData),
        colData = colData,
        design = ~ Age_group
      )
      dds <- DESeq(dds)
      
      res <- results(dds, contrast = c("Age_group", Age_group1, Age_group2))
      
      # 创建输出文件名
      base_name <- paste0(clean_filename(unique(rds$Region)), '_', clean_filename(f), '_', 
                         clean_filename(Age_group1), '_vs_', clean_filename(Age_group2))
      
      # 保存DESeq2结果
      write.csv(as.data.frame(res), 
                paste0(outdir, base_name, '_DESeq2_results.csv'), 
                row.names = TRUE)
      
      res <- na.omit(res)
      merge_TRADE <- TRADE(mode = "univariate", results1 = res)
      
      # 保存TRADE结果
      saveRDS(merge_TRADE, paste0(outdir1, base_name, '_trade_results.rds'))
      
      df1 <- as.data.frame(merge_TRADE$distribution_summary)
      df1 <- gather(df1)
      df1$Stage <- paste(Age_group1, Age_group2, sep = '_vs_')
      df1$Region <- unique(rds$Region)
      df1$subtype <- f
      
      write.csv(df1, paste0(outdir2, base_name, '_trade_results.csv'), 
                row.names = TRUE)
      
      message(paste("完成:", f, Age_group1, "vs", Age_group2))
      
    }, error = function(e) {
      message(paste("错误:", f, Age_group1, "vs", Age_group2, "-", e$message))
    })
  }
}