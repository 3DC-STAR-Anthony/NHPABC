# Load packages
library(ggplot2)
library(dplyr)
library(egg)
library(patchwork)
library(ggpubr)
library(ggbeeswarm)

###################################################################----Cell Proportion Age Group Boxplot (snRNA + snATAC side-by-side)----###################################################################
# Read input cell proportion table
df <- read.csv("./LTSR_cell_count/All_neuron_glia_prop_new_filter.csv", stringsAsFactors = FALSE)
df$Proportion <- df$n / df$Total_cell
# Filter low cell count samples
df <- df[df$Total_cell >= 100, ]

# Set age group order & color palette
age_order <- c('Young','Middle age','Old','Exceptionally old')
df$Age_group <- factor(df$Age_group, levels = age_order)
age_group_color <- c(
  'Young' = "#CC340B",
  'Middle age' = "#F08700",
  'Old' = "#9EC416",
  'Exceptionally old' = "#44C1EF"
)

# Plot parameters preset
point_size <- 0.4
box_alpha <- 0.5
jitter_width <- 0.2
linewidth_base <- 0.5
text_size <- 6
title_size <- 7

# Loop to plot for target cell type & region
celltype_l <- c("OPC")
region_l <- c("PM")
outdir <- "./Figure2_prop_boxplot/"

for (celltype in celltype_l) {
  # Filter data for target cell & region
  celltype_data <- df[df$Celltype == celltype & df$Region %in% region_l, ]
  
  # snRNA-seq boxplot
  p_rna <- ggplot(celltype_data[celltype_data$Batch == "snRNA-seq", ],
                   aes(x = Age_group, y = Proportion, fill = Age_group)) +
    geom_boxplot(alpha = box_alpha, outlier.shape = NA, linewidth = linewidth_base, position = position_dodge(0.8)) +
    geom_point(shape = 16, color = "black", size = point_size,
               position = position_jitterdodge(jitter.width = jitter_width, dodge.width = 0.8), alpha = 1) +
    scale_fill_manual(values = age_group_color) +
    labs(title = paste0(celltype, " Proportion (snRNA-seq)"), x = "", y = "Proportion", fill = "Age Group") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = linewidth_base),
      text = element_text(family = "sans", size = text_size, color = "black"),
      axis.text = element_text(family = "sans", size = text_size, color = "black"),
      axis.title = element_text(family = "sans", size = text_size, color = "black"),
      plot.title = element_text(family = "sans", size = title_size, color = "black"),
      legend.text = element_text(family = "sans", size = text_size, color = "black"),
      legend.title = element_text(family = "sans", size = text_size, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) + guides(color = "none")
  
  # snATAC-seq boxplot
  p_atac <- ggplot(celltype_data[celltype_data$Batch == "snATAC-seq", ],
                   aes(x = Age_group, y = Proportion, fill = Age_group)) +
    geom_boxplot(alpha = box_alpha, outlier.shape = NA, linewidth = linewidth_base, position = position_dodge(0.8)) +
    geom_point(shape = 16, color = "black", size = point_size,
               position = position_jitterdodge(jitter.width = jitter_width, dodge.width = 0.8), alpha = 1) +
    scale_fill_manual(values = age_group_color) +
    labs(title = paste0(celltype, " Proportion (snATAC-seq)"), x = "", y = "Proportion", fill = "Age Group") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = linewidth_base),
      text = element_text(family = "sans", size = text_size, color = "black"),
      axis.text = element_text(family = "sans", size = text_size, color = "black"),
      axis.title = element_text(family = "sans", size = text_size, color = "black"),
      plot.title = element_text(family = "sans", size = title_size, color = "black"),
      legend.text = element_text(family = "sans", size = text_size, color = "black"),
      legend.title = element_text(family = "sans", size = text_size, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) + guides(color = "none")
  
  # Unify y-axis range
  y_range <- range(celltype_data$Proportion, na.rm = TRUE)
  y_range <- c(y_range[1], y_range[2] * 1.1)
  p_rna <- p_rna + ylim(y_range)
  p_atac <- p_atac + ylim(y_range)
  
  # Combine two plots side by side
  combined_plot <- p_rna + p_atac + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")
  
  # Export individual PDF files
  file_rna <- paste0(outdir, region_l, "_", gsub("[^A-Za-z0-9-]", "_", celltype), "_proportion_comparison_snRNA.pdf")
  ggsave(file_rna, egg::set_panel_size(p_rna, width = unit(16, "mm"), height = unit(24, "mm")),
         width = 5, height = 5, units = "in", dpi = 300, limitsize = FALSE)
  
  file_atac <- paste0(outdir, region_l, "_", gsub("[^A-Za-z0-9-]", "_", celltype), "_proportion_comparison_snATAC.pdf")
  ggsave(file_atac, egg::set_panel_size(p_atac, width = unit(16, "mm"), height = unit(24, "mm")),
         width = 5, height = 5, units = "in", dpi = 300, limitsize = FALSE)
}


###################################################################----LR Interaction Number Boxplot----###################################################################
# Cell type color palette
cluster_color <- c() # Please add this information!

# Set factor levels for region and stage
region_order <- c()
stage_order <- c()
mat$region <- factor(mat$region, levels = region_order)
mat$stage <- factor(mat$stage, levels = stage_order)

# Create faceted boxplot
p_lr_box <- ggplot(mat, aes(x = region, y = Number, fill = bigcluster)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(aes(color = bigcluster), width = 0.2, size = 0.3) +
  scale_color_manual(values = cluster_color) +
  scale_fill_manual(values = cluster_color) +
  labs(x = "", y = "Number of interactions") +
  facet_wrap(~stage, ncol = 4) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6, family = "sans", color = "Black"),
    axis.text.y = element_text(size = 6, family = "sans", color = "Black"),
    axis.title.y = element_text(size = 6, family = "sans", color = "Black"),
    legend.title = element_text(size = 6, color = "Black"),
    legend.text = element_text(size = 6, color = "Black"),
    legend.position = "left",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

###################################################################----DARE Percentage Stage Comparison Boxplot----###################################################################
# Read input data
sdare_df <- read.csv("./sdare_percent_table.csv", stringsAsFactors = FALSE)

# Stage color palette
stage_col <- c() # Please fill your stage color mapping here

# Plot boxplot + quasirandom scatter
p_sdare_box <- ggplot(sdare_df, aes(x = stage, y = percent, color = stage)) +
  geom_quasirandom(dodge.width = 0.8, size = 0.02) +
  geom_boxplot(outlier.shape = NA, width = 0.8, fill = NA) +
  scale_color_manual(values = stage_col) +
  ylim(c(0, 110)) +
  xlab("") +
  ylab("Percentage") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1, family = "sans", color = "black"),
    axis.text.y = element_text(size = 6, family = "sans", color = "black"),
    axis.title.x = element_text(size = 6, family = "sans", color = "black"),
    axis.title.y = element_text(size = 6, family = "sans", color = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 6, hjust = 0.5, family = "sans", color = "black")
  )

# Export figure
ggsave(
  "NHPABC_sDAREs_subtype_percentage_boxplot.pdf",
  egg::set_panel_size(p_sdare_box, width = unit(21, "mm"), height = unit(16, "mm")),
  height = 4, width = 15, limitsize = FALSE
)

###################################################################----DARE Up/Down Trend Percentage Boxplot----###################################################################
# Count table for each stage
df_pdare <- as.data.frame(table(pdare_all$label, pdare_all$coef_trend), stringsAsFactors = FALSE)
df_early <- as.data.frame(table(early_all$label, early_all$trend), stringsAsFactors = FALSE)
df_late <- as.data.frame(table(late_all$label, late_all$trend), stringsAsFactors = FALSE)
df_verylate <- as.data.frame(table(verylate_all$label, verylate_all$trend), stringsAsFactors = FALSE)

# Mark stage info
df_pdare$stage <- "pDAREs"
df_early$stage <- "Early"
df_late$stage <- "Late"
df_verylate$stage <- "Very late"

# Merge all data frames
df_dare_list <- list(df_pdare, df_early, df_late, df_verylate)
df_dare <- Reduce(rbind, df_dare_list)

# Calculate percentage within each cell subtype
df_dare <- df_dare %>%
  group_by(Var1) %>%
  mutate(
    var_total = sum(Freq),
    percent = Freq / var_total * 100
  ) %>%
  ungroup()

# Set stage display order
df_dare$stage <- factor(df_dare$stage, levels = c("pDAREs","Early","Late","Very late"))

# Trend color palette
trend_color <- c() # Please fill your Up/Down color mapping here

# Draw boxplot + scattered points
p_dare_trend_box <- ggplot(df_dare, aes(x = stage, y = percent, color = Var2)) +
  geom_quasirandom(dodge.width = 0.8, size = 0.02) +
  geom_boxplot(outlier.shape = NA, width = 0.8, fill = NA) +
  scale_color_manual(values = trend_color, name = "Trend") +
  ylim(c(0, 120)) +
  xlab("") +
  ylab("Percentage") +
  ggtitle("DAREs") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1, family = "sans", color = "black"),
    axis.text.y = element_text(size = 6, family = "sans", color = "black"),
    axis.title.x = element_text(size = 6, family = "sans", color = "black"),
    axis.title.y = element_text(size = 6, family = "sans", color = "black"),
    legend.title = element_text(size = 6, family = "sans", color = "black"),
    legend.text = element_text(size = 6, family = "sans", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 8, hjust = 0.5, family = "sans", color = "black")
  )

# Export PDF figure
ggsave(
  "NHPABC_DAREs_trend_percentage_boxplot.pdf",
  egg::set_panel_size(p_dare_trend_box, width = unit(38, "mm"), height = unit(16, "mm")),
  height = 4, width = 15, limitsize = FALSE
)

###################################################################----LR Ligand & Receptor Enrichment Score Boxplot----###################################################################
# Read input CellChat LR table
charlr <- fread("./NHPABC_CellChat_truncatedMean_final_ARLR.csv")

# Step1: Get all unique cell subtypes from source & target labels
lc_col_an <- data.frame(celltype = unique(c(charlr$source_label, charl r$target_label)))
rownames(lc_col_an) <- lc_col_an$celltype

# Step2: Build ligand count matrix (source cell type count per pathway)
lg_mat <- matrix(0, nrow = length(unique(charlr$pathway_name)), ncol = nrow(lc_col_an))
colnames(lg_mat) <- rownames(lc_col_an)
rownames(lg_mat) <- unique(charlr$pathway_name)
for(f in rownames(lg_mat)){
  df_sub <- charl r[charlr$pathway_name == f, ]
  for(i in colnames(lg_mat)){
    cnt <- nrow(df_sub[df_sub$source_label == i, ])
    lg_mat[f, i] <- cnt
  }
}
# Convert raw count to proportion matrix (row sum = 1)
lg_mat_per_subtype <- lg_mat
for(i in colnames(lg_mat_per_subtype)){
  lg_mat_per_subtype[,i] <- lg_mat_per_subtype[,i] / rowSums(lg_mat)
}

# Step3: Build receptor count matrix (target cell type count per pathway)
rp_mat <- matrix(0, nrow = length(unique(charlr$pathway_name)), ncol = nrow(lc_col_an))
colnames(rp_mat) <- rownames(lc_col_an)
rownames(rp_mat) <- unique(charlr$pathway_name)
for(f in rownames(rp_mat)){
  df_sub <- charl r[charlr$pathway_name == f, ]
  for(i in colnames(rp_mat)){
    cnt <- nrow(df_sub[df_sub$target_label == i, ])
    rp_mat[f, i] <- cnt
  }
}
# Convert raw count to proportion matrix
rp_mat_per_subtype <- rp_mat
for(i in colnames(rp_mat_per_subtype)){
  rp_mat_per_subtype[,i] <- rp_mat_per_subtype[,i] / rowSums(rp_mat)
}

# Step4: Reshape matrix to long data frame for plotting
# Ligand table
dfl <- as.data.frame(t(lg_mat_per_subtype[unique(charlr$pathway_name), ]))
dfl$subtype <- rownames(dfl)
dfl <- gather(dfl, key = pathway, value = ligand, -subtype)
# Receptor table
dfr <- as.data.frame(t(rp_mat_per_subtype[unique(charlr$pathway_name), ]))
dfr$subtype <- rownames(dfr)
dfr <- gather(dfr, key = pathway, value = receptor, -subtype)
# Merge ligand & receptor values
dfmerge <- dfl
dfmerge$receptor <- dfr$receptor

# Step5: Reusable boxplot drawing function
single_barplot <- function(f, dt, savepath, width = 18, height = 15){
  df_sub <- dt[dt$pathway == f, ]
  df_sub$bigcluster <- factor(df_sub$bigcluster, levels = c('Ex','Inh','Other neuron','Ast','Mic','ODC','OPC','VS','Ependymal'))
  max_lim <- max(c(df_sub$ligand, df_sub$receptor)) + 0.01
  
  # Ligand enrichment boxplot
  p_lig <- ggplot(df_sub, aes(x = bigcluster, y = ligand, fill = bigcluster)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = region), size = 0.05, width = 0.25) +
    scale_fill_manual(values = ano_col$Bigcluster, name = "Cell type") +
    scale_color_manual(values = ano_col$Region, name = "Region") +
    theme_bw() +
    xlab("") +
    ylab("Ligand enrichment score") +
    ggtitle(paste0(f, " LR pathway enrichment")) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max_lim)) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 6, family = "sans", color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6, family = "sans", color = "black"),
      legend.text = element_text(size = 6, family = "sans", color = "black"),
      legend.title = element_text(size = 6, family = "sans", color = "black"),
      legend.position = "right",
      legend.box = "horizontal",
      legend.spacing.x = unit(0.5, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5, family = "sans", color = "black")
    )
  ggsave(paste0(savepath, f, "_L_enrich_boxplot.pdf"),
         egg::set_panel_size(p_lig, width = unit(width, "mm"), height = unit(height, "mm")),
         height = 5, width = 5, limitsize = FALSE)
  
  # Receptor enrichment boxplot
  p_rec <- ggplot(df_sub, aes(x = bigcluster, y = receptor, fill = bigcluster)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = region), size = 0.05, width = 0.25) +
    scale_fill_manual(values = ano_col$Bigcluster, name = "Cell type") +
    scale_color_manual(values = ano_col$Region, name = "Region") +
    theme_bw() +
    xlab("") +
    ylab("Receptor enrichment score") +
    ggtitle(paste0(f, " receptor pathway enrichment")) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max_lim)) +
    theme(
      axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, family = "sans", color = "black"),
      axis.text.y = element_text(size = 6, family = "sans", color = "black"),
      axis.title.x = element_text(size = 6, family = "sans", color = "black"),
      axis.title.y = element_text(size = 6, family = "sans", color = "black"),
      legend.text = element_text(size = 6, family = "sans", color = "black"),
      legend.title = element_text(size = 6, family = "sans", color = "black"),
      legend.position = "right",
      legend.box = "horizontal",
      legend.spacing.x = unit(0.5, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank()
    )
  ggsave(paste0(savepath, f, "_R_enrich_boxplot.pdf"),
         egg::set_panel_size(p_rec, width = unit(width, "mm"), height = unit(height, "mm")),
         height = 5, width = 5, limitsize = FALSE)
}

# Define output path & target pathway list
savepath <- "./Figure5/Single_LR_enrich_box/"
target_pathways <- c("EPHA", "NRXN", "EGF", "TGFB", "PDGF")
# Loop to generate boxplot for each pathway
for(pathway_name in target_pathways){
  single_barplot(f = pathway_name, dt = dfmerge, savepath = savepath, width = 18, height = 15)
}
