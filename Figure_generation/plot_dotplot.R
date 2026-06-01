# Load packages
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(egg)

###################################################################----RNA Marker Gene Dotplot----###################################################################
# Set fixed cell type order
Cell_level <- c()

# Manually curated global RNA marker gene list
global_rna_marker <- c()

# Load Seurat object and set cell type factor level
scrna_obj <- readRDS("./scrna_seurat.rds")
scrna_obj@meta.data$cell_type <- factor(scrna_obj@meta.data$cell_type, levels = Cell_level)

# Generate marker dotplot
p_marker_dot <- DotPlot(
  scrna_obj,
  group.by = 'cell_type',
  features = global_rna_marker,
  scale = FALSE,
  dot.scale = 3,
  scale.by = "radius",
  dot.min = -3
) +
scale_color_gradientn(values = seq(0,1,0.01), colours = brewer.pal(9,"Reds")) +
theme_bw() +
theme(
  axis.text.x = element_text(family = "sans", angle = 45, size = 6, color = "black", hjust = 1, vjust = 1),
  axis.text.y = element_text(family = "sans", size = 6, color = "black"),
  legend.title = element_text(family = "sans", size = 6, color = "black"),
  legend.text = element_text(family = "sans", size = 6, color = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
) +
ylab("") +
xlab("") +
scale_y_discrete(limits = rev(Cell_level))

# Save PDF figure
ggsave(
  "NHPABC_RNA_marker_dotplot.pdf",
  egg::set_panel_size(p_marker_dot, width = unit(120, "mm"), height = unit(30, "mm")),
  width = 15, height = 3, limitsize = FALSE
)

###################################################################----LTSR Fold Change Dotplot----###################################################################
# Read input data
merged_data <- read.csv("merged_data_LTSR_prop.csv", stringsAsFactors = FALSE)

# Filter significant records and export table
merged_data_filter <- merged_data[merged_data$LFSR <= 0.01 & abs(merged_data$Estimate) >= 0.005, ]

# Prepare plotting dataset
all_ltsr <- merged_data
all_ltsr$X <- all_ltsr$Celltype
all_ltsr$Age_ltsr <- all_ltsr$LFSR
all_ltsr$Age_postmean <- all_ltsr$Estimate

# Fixed plotting order vectors
region_lev <- c()
celltype_lev <- c()

# Predefine size breakpoints and color gradient limits
size_breaks <- c(">0.9999",">0.99",">0.9","<=0.9")
size_values <- c(3,2.5,1,0.5)
color_midpoint <- 0
range_num <- c(-0.9,0.9)

# Loop plot by each age stage
stages <- unique(all_ltsr$Stage)
plot_list <- list()

for (stage in stages) {
  stage_data <- all_ltsr %>% filter(Stage == stage) %>% mutate(
    size_level = case_when(
      Age_ltsr > 0.9999 ~ ">0.9999",
      Age_ltsr > 0.99 ~ ">0.99",
      Age_ltsr > 0.9 ~ ">0.9",
      TRUE ~ "<=0.9"
    ),
    size_level = factor(size_level, levels = size_breaks)
  )
  # Set factor order
  stage_data$Region <- factor(stage_data$Region, levels = region_lev)
  stage_data$X <- factor(stage_data$X, levels = celltype_lev)
  
  # Draw dot plot
  p <- ggplot(stage_data, aes(x = Region, y = X, stroke = 0)) +
    geom_point(aes(size = size_level, color = Age_postmean), alpha = 0.7) +
    scale_color_gradient2(
      low = "#1160A0", mid = "lightgrey", high = "red",
      midpoint = color_midpoint, limit = range_num, name = "Fold Change"
    ) +
    scale_size_manual(values = size_values, breaks = size_breaks, name = "Age_LTSR") +
    labs(x = "Celltype", y = "Region", size = "ltsr", title = stage) +
    theme_bw() +
    theme(
      text = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    coord_flip()
  plot_list[[stage]] <- p
}

# Export single stage example figure (p_prop stage)
p_prop <- plot_list[["p_prop"]]
ggsave(
  "Neuron_proportion_LTSR.pdf",
  egg::set_panel_size(p_prop, width = unit(160, "mm"), height = unit(24, "mm")),
  width = 10, height = 10, units = "in", dpi = 300, limitsize = FALSE
)


###################################################################----Ligand/Receptor Age Coefficient Bubble Dotplot----###################################################################
# Reusable bubble dotplot function
create_expression_bubbleplot <- function(dfmt, f, text_size = 6, type = "ligand") {
  # Dynamic y axis position for ligand/receptor
  y_position <- ifelse(type == "receptor", "left", "right")
  
  p <- ggplot(dfmt, aes(x = label, y = gene)) +
    geom_point(aes(size = pct.exp, fill = coef, color = Sig), shape = 21, stroke = 0.5) +
    scale_size("Percentage expression", range = c(1, 2)) +
    scale_fill_gradient2(low = "royalblue", mid = "white", high = "#FF4040", name = "Age coefficient") +
    scale_color_manual(values = c("Non-sig" = "grey75", "Sig" = "BLACK"), name = "Significant") +
    scale_y_discrete(position = y_position) +
    theme_cowplot() +
    xlab("") +
    ylab("") +
    theme_bw() +
    ggtitle(f) +
    theme(
      axis.text.x = element_text(family = "sans", size = text_size, angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
      axis.text.y = element_text(family = "sans", size = text_size, colour = "black"),
      legend.text = element_text(family = "sans", colour = "black", size = text_size),
      legend.title = element_text(family = "sans", colour = "black", size = text_size),
      legend.key.height = unit(10, "mm"),
      legend.key.width = unit(2, "mm"),
      plot.title = element_text(family = "sans", size = text_size + 2, hjust = 0.5),
      panel.grid.major = element_blank()
    )
  return(p)
}

# Load & filter input ligand-receptor data
lr_dot_exp <- read.csv("./LR_coef_expression_table.csv", stringsAsFactors = FALSE)
# Predefined selected cell type labels
label_sel <- c()

# Subset ligand data
gene.ligand <- unique(lr_dot_exp$gene[lr_dot_exp$gene.ligand])
dfmt <- lr_dot_exp[lr_dot_exp$gene %in% gene.ligand & lr_dot_exp$label %in% label_sel, ]
# Set fixed factor order
dfmt$label <- factor(dfmt$label, label_sel)
dfmt$gene <- factor(dfmt$gene, rev(gene.ligand))

# Generate ligand bubble dotplot
p2 <- create_expression_bubbleplot(dfmt, f = "ligand", text_size = 6, type = "ligand")

# Export PDF figure
ggsave(
  paste0("XXX", "_ligand", "_dotplot.pdf"),
  egg::set_panel_size(p2, width = unit(50, "mm"), height = unit(15, "mm")),
  height = 5, width = 5, limitsize = FALSE
)
