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

# Load Seurat object
rds <- readRDS('~/scRNA/Each_single_region/')

# Load protein coding gene list
protein_cd <- read.csv('~/scRNA/protein_coding_gene_list.csv')

# Filter out LOC1 genes
anno_gene <- unique(protein_cd$Gene)[grep('LOC1', unique(protein_cd$Gene), invert = TRUE)]
filter <- subset(rds, features = anno_gene)

# Load DEG list
deg_list <- read.csv('~/Figure3_DEG/Revise_DEG_LIST/Revised_DEG_all_filter.csv')

# Define output directories
outdir <- '~/Figure3_DEG/Trade2/Trade_DEseq2/'
outdir1 <- '~/Figure3_DEG/Trade2/Trade_result_list/'
outdir2 <- '~/Figure3_DEG/Trade2/Trade_result_table/'

# Get unique regions
Region_list <- unique(rds$Region)
subtypes <- unique(deg_list$subtype[deg_list$Region == Region_list])

# Define age group combinations for comparison
age_combinations <- list(
  c('Young', 'Middle age'),
  c('Middle age', 'Old'),
  c('Old', 'Exceptionally old')
)

# Function: Clean filename by removing invalid characters
clean_filename <- function(name) {
  # Replace special characters
  clean_name <- gsub('[+/\\\\:*?"<>| ]', '_', name)
  # Remove leading/trailing underscores
  clean_name <- gsub('^_+|_+$', '', clean_name)
  # Replace multiple consecutive underscores with single
  clean_name <- gsub('_+', '_', clean_name)
  return(clean_name)
}

# Process each subtype and age combination
for(f in subtypes) {
  for(age_pair in age_combinations) {
    Age_group1 <- age_pair[1]
    Age_group2 <- age_pair[2]
    
    # Subset data for current subtype
    filter1 <- subset(filter, subtype == f)
    filter2 <- subset(filter1, Age_group == Age_group1 | Age_group == Age_group2)
    
    # Check if there is enough data
    if(nrow(filter2@meta.data) == 0) {
      message(paste("Skipping:", f, Age_group1, "vs", Age_group2, "- no data"))
      next
    }
    
    # Aggregate expression by Age_group and Sample
    aggr <- AggregateExpression(filter2, 
                                assays = 'RNA',
                                group.by = c("Age_group", "Sample"),
                                return.seurat = FALSE)
    
    # Check if there are enough age groups for DESeq2 analysis
    if(length(unique(filter2$Age_group)) < 2) {
      message(paste("Skipping:", f, Age_group1, "vs", Age_group2, "- insufficient age groups"))
      next
    }
    
    av <- as.data.frame(aggr[[1]])
    countData <- av
    
    # Check if there are enough columns
    if(ncol(countData) < 2) {
      message(paste("Skipping:", f, Age_group1, "vs", Age_group2, "- insufficient samples"))
      next
    }
    
    # Prepare column data
    colData <- str_split_fixed(colnames(countData), '_', n = 2) %>% 
      as.data.frame() %>% 
      setNames(c('Age_group', 'Sample'))
    colData$id <- colnames(countData)
    rownames(colData) <- colData$id
    
    # DESeq2 analysis
    tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = round(countData),
        colData = colData,
        design = ~ Age_group
      )
      dds <- DESeq(dds)
      
      # Get differential expression results
      res <- results(dds, contrast = c("Age_group", Age_group1, Age_group2))
      
      # Create output file name
      base_name <- paste0(clean_filename(unique(rds$Region)), '_', clean_filename(f), '_', 
                         clean_filename(Age_group1), '_vs_', clean_filename(Age_group2))
      
      # Save DESeq2 results
      write.csv(as.data.frame(res), 
                paste0(outdir, base_name, '_DESeq2_results.csv'), 
                row.names = TRUE)
      
      # Remove NA values
      res <- na.omit(res)
      
      # Perform TRADE analysis
      merge_TRADE <- TRADE(mode = "univariate", results1 = res)
      
      # Save TRADE results
      saveRDS(merge_TRADE, paste0(outdir1, base_name, '_trade_results.rds'))
      
      # Process TRADE results
      df1 <- as.data.frame(merge_TRADE$distribution_summary)
      df1 <- gather(df1)
      df1$Stage <- paste(Age_group1, Age_group2, sep = '_vs_')
      df1$Region <- unique(rds$Region)
      df1$subtype <- f
      
      # Save TRADE summary
      write.csv(df1, paste0(outdir2, base_name, '_trade_results.csv'), 
                row.names = TRUE)
      
      message(paste("Completed:", f, Age_group1, "vs", Age_group2))
      
    }, error = function(e) {
      message(paste("Error:", f, Age_group1, "vs", Age_group2, "-", e$message))
    })
  }
}
