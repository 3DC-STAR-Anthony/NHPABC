# Load required packages
library(ArchR)
library(egg)

###################################################################----Peak-to-Gene Links Side-by-Side Heatmap (scATAC + scRNA)----###################################################################
# Load ArchR project
projHeme5 <- loadArchRProject(path = "./ArchRProject/")

# Generate paired ATAC-RNA peak-to-gene heatmap
p_peak2gene <- plotPeak2GeneHeatmap(
  ArchRProj = projHeme5,
  groupBy = "Clusters2"
)

# Save output PDF
ggsave(
  "NHPABC_Peak2Gene_dual_heatmap.pdf",
  egg::set_panel_size(p_peak2gene, width = unit(100, "mm"), height = unit(80, "mm")),
  width = 12, height = 10, limitsize = FALSE
)

# Load packages
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(egg)

###################################################################----DEG Module Overlap Odds Ratio Heatmap----###################################################################
# Input data
df1 <- read.csv("./DGE_table.csv", stringsAsFactors = FALSE)
protein_coding_gene_list <- read.csv("protein_coding_gene_list.csv", stringsAsFactors = F)

# Filter DGE table
df1 <- df1[abs(df1$abs) >= 0.005, ]
df1$Gene_2 <- with(df1, paste0(Gene, "|", Trend))

# Initialize overlap storage
overlap_df <- df1
overlap_df$slspldap <- NULL
overlap_merge <- NULL
gene <- protein_coding_gene_list
modules <- df1
mods <- unique(modules$label)
genome.size <- length(unique(gene$Gene))

# Loop to calculate gene overlap statistics
for (f in unique(df1$label)) {
  over_gene <- subset(df1, label == f)$Gene_2
  cur_overlap_df <- do.call(rbind, lapply(mods, function(cur_mod) {
    cur_genes <- modules %>% subset(label == cur_mod) %>% .$Gene_2
    cur_overlap <- testGeneOverlap(newGeneOverlap(over_gene, cur_genes, genomeSize = genome.size))
    data.frame(
      'odds.ratio' = cur_overlap@odds.ratio,
      'pval' = cur_overlap@pval,
      'Jaccard' = cur_overlap@Jaccard,
      'size_intersection' = length(cur_overlap@intersection),
      'module' = cur_mod
    )
  }))
  cur_overlap_df$fdr <- p.adjust(cur_overlap_df$pval, method = 'fdr')
  cur_overlap_df$shape <- ifelse(cur_overlap_df$fdr <= 0.05, 4, 21)
  cur_overlap_df$module <- factor(as.character(cur_overlap_df$module), levels = as.character(mods))
  cur_overlap_df$label <- f
  cur_overlap_df$Subtype <- unique(subset(df1, label == f)$Subtype)
  overlap_merge <- rbind(overlap_merge, cur_overlap_df)
}

# Process odds ratio extreme values
max_odds <- max(subset(overlap_merge, odds.ratio != Inf)$odds.ratio)
min_odds <- min(subset(overlap_merge, odds.ratio != 0)$odds.ratio)
overlap_merge$odds.ratio[overlap_merge$odds.ratio == Inf] <- max_odds
overlap_merge$odds.ratio[overlap_merge$odds.ratio == 0] <- min_odds

# Mark significance
overlap_merge$sig <- ifelse(overlap_merge$fdr > 0.01, 'ns', 'sig')
loss <- subset(overlap_merge, fdr > 0.01)

# Log transform odds ratio
overlap_merge$log10OR <- log10(overlap_merge$odds.ratio)

# Cast long table to wide matrix for heatmap
heat1 <- overlap_merge %>% dplyr::select(c('log10OR', 'module', 'label')) %>% spread(module, log10OR)
rownames(heat1) <- heat1$label
heat1 <- dplyr::select(heat1, -label)

# Generate pheatmap
bk <- c(seq(min(heat1, na.rm = TRUE), 0, by = 0.01), seq(0, max(heat1, na.rm = TRUE), by = 0.01))
p1 <- pheatmap(
  heat1,
  na_col = 'grey',
  clustering_method = 'ward.D2',
  color = c(
    colorRampPalette(c("#9A3F8D","#3F4697"))(50),
    colorRampPalette(c("#3F4697","#1C142F"))(100),
    colorRampPalette(c("#1C142F","#FAE577"))(150)
  ),
  border = NA,
  number_color = "white",
  cellwidth = 5,
  cellheight = 5,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 5,
  fontsize_col = 5,
  angle_col = "90",
  file = "NHPABC_DEG_module_overlap_oddsratio_heatmap.pdf"
)

###################################################################----Pathway DEG Coefficient Heatmap----###################################################################
# Custom matrix fill function
fillARGsval <- function(mt, arg, celltype, gene, fil = TRUE){
    for(f in celltype){
        for(i in gene){
            a <- arg[arg$label == f & arg$gene == i, "coef"]
            if(dim(arg[arg$label == f & arg$gene == i,])[1] == 0){
                mt[i, f] <- NA
            }else{
                mt[i, f] <- a
            }
        }
    }
    if(fil == TRUE){
        mt <- mt[which(rowSums(is.na(mt)) < length(celltype)), ]
    }
    return(mt)
}

# Load input gene coefficient table
df_gene_all <- read.csv("./Selected_gene.csv", stringsAsFactors = FALSE)
lev$label <- unique(df_gene_all$label)

# Initialize empty matrix
mt <- matrix(NA, ncol = length(lev$label), nrow = length(rownames(df_gene_all)))
colnames(mt) <- lev$label
rownames(mt) <- df_gene_all$gene

# Fill matrix with gene coef values
up <- fillARGsval(mt, df_gene_all, lev$label, df_gene_all$gene, fil = FALSE)
up[is.na(up)] <- 0
colnames(up) <- lev$label

# Prepare row & column annotation data
anno_row <- df_gene_all[, c(2, 1)]
rownames(anno_row) <- df_gene_all$gene
anno_col <- lev[lev$label %in% colnames(up), ]
rownames(anno_col) <- anno_col$label
anno_col <- anno_col[, c(-1, -4)]
colnames(anno_col) <- c("Region", "Bigcluster")

# Extract subtype labels for row labels
lev$subtype <- str_extract(lev$label, ".*(?=[_]$)")

# Color breakpoints setup
bk <- c(seq(-max(abs(up)), -0.0001, by = 0.0001), seq(0, max(abs(up)), by = 0.0001))

# Generate pheatmap
p_coef_heat <- pheatmap(
    up,
    main = paste0("Ligand and Receptor"),
    breaks = bk,
    border_color = NA,
    fontsize = 6,
    cellwidth = 3.2,
    cellheight = 7.8,
    annotation_col = anno_col,
    annotation_row = anno_row,
    annotation_colors = anno_col,
    cluster_rows = F,
    cluster_cols = F,
    labels_row = lev$subtype,
    #gaps_row = c(6,8,10,16,18,22,25,34,43,48,58,66,72),
    #gaps_col = c(6,8,10,12,16,19,28,37,42,52,56,64),
    angle_col = 90,
    color = c(
        colorRampPalette(c("blue", "white"))(length(bk)/2),
        colorRampPalette(c("white", "red"))(length(bk)/2)
    ),
    filename = "NHPABC_DEG_coef_heatmap.pdf"
)

###################################################################----cCRE Peak CPM Heatmap----###################################################################
# Load input data
cep <- read.csv("./cCRE_peak_subtype.csv", stringsAsFactors = FALSE)
allcpm <- read.csv("./cCRE_cpm_matrix.csv", stringsAsFactors = FALSE)
lev <- read.csv("./celltype_anno_info.csv", stringsAsFactors = FALSE)

# Deduplicate peaks
cep4_un <- cep4 %>% distinct(peaks, .keep_all = TRUE)
allcpm_un <- allcpm %>% distinct(peaks, .keep_all = TRUE)

# Initialize empty matrix
mat <- matrix(NA, nrow = length(lev$label), ncol = length(cep_un$peaks))
rownames(mat) <- lev$label
colnames(mat) <- cep4_un$peaks

# Fill matrix with CPM values
for(f in rownames(mat)){
  ccre_raw <- allcpm_un[allcpm_un$label == f, ]
  mat[f, ] <- ccre_raw$cpm[match(colnames(mat), rownames(ccre_raw))]
}

# Replace zero with NA and log10 transform
mat[mat == 0] <- NA
matlog <- log10(mat)

# Prepare row annotation and color list
anno_row <- lev[, c("region", "Bigcluster")]
rownames(anno_row) <- lev$label
anno_color <- list(
  Region = c(
    PFC = '#1B9E77',
    Hf = '#D95F02',
    Str = '#7570B3',
    Tha = '#E7298A',
    Hyt = '#66A61E',
    MB = '#E6AB02',
    PM = '#A6761D',
    Cer = '#666666'
  ),
  Type = c(
    'Ex' = "#8FE0C0",
    'Inh' = "#4E4179",
    'other neuron' = "#8FA072",
    'ODC' = "#55978A",
    'OPC' = "#3E3E39",
    'Ast' = "#A73E1B",
    'Mic' = "#49499C",
    'Ependymal' = "#B57014",
    'VS' = "#808080"
  )
)

# Plot cCRE heatmap
p_ccre_heat <- pheatmap(
  matlog,
  color = c(
    colorRampPalette(c("#441753","#C53C88"))(80),
    colorRampPalette(c("#C53C88","#D6DF28"))(30),
    colorRampPalette(c("#D6DF28","#FFED1F"))(40)
  ),
  scale = "none",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  annotation_row = anno_row,
  annotation_colors = anno_color,
  angle_col = 90,
  fontsize = 15,
  cellwidth = 30,
  legend = FALSE,
  filename = "NHPABC_cCRE_heatmap.pdf"
)

# Load packages
library(dplyr)
library(pheatmap)
library(stringr)

###################################################################----cCRE Motif Enrichment Heatmap (Common + Specific Motifs)----###################################################################
# Custom functions
create_enrichment_matrix <- function(df, celltype, feature, trend = "NULL"){
  result <- matrix(0, nrow = length(unique(celltype)), ncol = length(unique(feature)), dimnames = list(unique(celltype), unique(feature)))
  idx <- cbind(match(df$celltype, rownames(result)), match(df$feature, colnames(result)))
  result[idx] <- df$log10Padj
  if(trend == "Down") result <- -result
  return(result)
}
find_common_motifs <- function(mat, threshold_rows = 10, p_value = 0.05){
  threshold_value <- -log10(p_value)
  motif_stats <- apply(mat, 2, function(col) sum(col >= threshold_value, na.rm = TRUE))
  return(names(motif_stats)[motif_stats >= threshold_rows])
}

# Load data
jasper2024_list <- read.csv("./jasper2024_motif_gene_list.csv", stringsAsFactors = FALSE)

# Build matrix
mt_ccre <- create_enrichment_matrix(jasper2024_list, lev$label, jasper2024_list$feature)
mt_ccre_filter <- mt_ccre[, colSums(mt_ccre != 0) > 0]

# Get common/specific motifs
motif_common_list <- find_common_motifs(mt_ccre_filter, 30, 0.05)
mt_ccre_common <- mt_ccre_filter[, motif_common_list]
mt_ccre_specific <- mt_ccre_filter[, !colnames(mt_ccre_filter) %in% motif_common_list]

# Sort
mt_ccre_common <- mt_ccre_common[, names(apply(mt_ccre_common, 2, which.max)[order(apply(mt_ccre_common, 2, which.max))])]
mt_ccre_specific <- mt_ccre_specific[, names(apply(mt_ccre_specific, 2, which.max)[order(apply(mt_ccre_specific, 2, which.max))])]

# Save and reload motif order
motif_comb <- c(motif_common_list, colnames(mt_ccre_specific))
write.csv(data.frame(name = motif_comb), "cCRE_motif_heatmap_level.csv", row.names = FALSE)
motif_lev_new <- read.csv("cCRE_motif_heatmap_level.csv", stringsAsFactors = FALSE)$name
data_to_plot <- mt_ccre_filter[, motif_lev_new]
colnames(data_to_plot)[!colnames(data_to_plot) %in% colnames(mt_ccre_specific)] <- ""

# Plot heatmap
bk <- seq(0, 100, by = 0.01)
p_final_motif <- pheatmap(
  data_to_plot,
  breaks = bk,
  fontsize = 6,
  border_color = "white",
  cellwidth = 0.74,
  cellheight = 0.82,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  angle_col = 90,
  cut_cols = 4,
  show_colnames = TRUE,
  color = c(colorRampPalette(c("#0B0405FF","#40498EFF"))(length(bk)/5),
            colorRampPalette(c("#40498EFF","#38AAACFF"))(2*length(bk)/5),
            colorRampPalette(c("#38AAACFF","#DEF5E5FF"))(2*length(bk)/5)),
  filename = "All_cCRE_Common_Specific_motif_for_figure.pdf"
)

# Load required packages
library(pheatmap)
library(patchwork)
library(dplyr)

###################################################################----GWAS cCRE Conservation & Non-Conservation Select Trait Heatmap----###################################################################
# Load input matrix data (Here we use the calibrated p-value as the numerical value to draw the heatmap.)
heat_con <- readRDS("./cCRE_conservation_trait_heatmap_matrix.rds")
heat_noncon <- readRDS("./cCRE_non_conservation_trait_heatmap_matrix.rds")
gwas_lev <- read.csv("./GWAS_trait_accession.csv", stringsAsFactors = FALSE)
rownames(gwas_lev) <- gwas_lev$Study_accession

# Trait grouping list
trait_lev <- c(
  # Neurological disorders
  "Neurological disorders",
  "Neurological features",
  "Other"
)
# Manually define GWAS trait grouping table
df_ano <- data.frame(
  GWAS = c(rep("Neurological disorders",4), rep("Neurological features",6), rep("Other",3)),
  time = seq_len(nrow(gwas_lev))
)
rownames(df_ano) <- gwas_lev$traitlev

# Subset matrix to selected GWAS traits
hm <- heat_con[rownames(df_ano), ]
hms <- heat_noncon[rownames(df_ano), ]

# Row annotation color config
annotation_col <- list(
  Bigcluster = c(
    Ex = '#8FE0C0', Inh = '#4E4179', `Other neuron` = '#FA8072',
    Ast = '#A73E1B', Mic = '#49499C', ODC = '#55978A',
    OPC = '#3E3E39', VS = '#808080', Ependymal = '#B57014'
  ),
  Region = c(
    PFC = '#1B9E77', Hf = '#D95F02', Str = '#7570B3', Tha = '#E7298A',
    Hyt = '#66A61E', MB = '#E6AB02', PM = '#A6761D', Cer = '#666666'
  )
)

# Color breakpoints
bk <- seq(0, max(abs(hm)), by = 0.0001)

# Plot Conservation GWAS heatmap
p1 <- pheatmap(
  hm,
  display_numbers = matrix(ifelse(hm > -log10(0.05), "", ""), nrow = nrow(hm)),
  number_color = "black",
  main = "GWAS (Conservation)",
  breaks = bk,
  fontsize = 6,
  border_color = NA,
  annotation_col = lev[, c(2,3)],
  annotation_row = df_ano,
  annotation_colors = annotation_col,
  cluster_cols = F,
  cluster_rows = F,
  gaps_row = c(4),
  gaps_col = c(20,43,51,59,67,75),
  cellwidth = 4.7,
  cellheight = 6,
  angle_col = 90,
  color = colorRampPalette(c("white","red"))(length(bk)),
  filename = "cCRE_CONSERVATION_GWAS_select_clustering_by_manual.pdf"
)

# Plot Non-Conservation GWAS heatmap
p2 <- pheatmap(
  hms,
  display_numbers = matrix(ifelse(hms > -log10(0.05), "", ""), nrow = nrow(hms)),
  number_color = "black",
  main = "GWAS (Non-Conservation)",
  breaks = bk,
  fontsize = 6,
  border_color = NA,
  annotation_col = lev[, c(2,3)],
  annotation_row = df_ano,
  annotation_colors = annotation_col,
  cluster_cols = F,
  cluster_rows = F,
  gaps_row = c(4),
  gaps_col = c(20,43,51,59,67,75),
  cellwidth = 4.7,
  cellheight = 6,
  angle_col = 90,
  color = colorRampPalette(c("white","red"))(length(bk)),
  filename = "cCRE_Non_Conservation_GWAS_select_clustering_by_manual.pdf"
)

# Combine two heatmaps side-by-side and export full figure
p <- plot_grid(p1$gtable, p2$gtable, ncol = 2, align = "hv", rel_heights = 1, rel_widths = 1, labels = c("Conservation","Non-Conservation"), label_size = 20, label_x = 0)
ggsave(
  "All_cCRE_conser_GWAS_select_clustering_by_manual.pdf",
  p,
  width = 15, height = 8, limitsize = FALSE
)                
