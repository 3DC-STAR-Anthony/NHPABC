library(ggplot2)
library(cowplot)
library(egg)
library(Seurat)
library(ArchR)

# Load input data (modify path as needed)
input_file <- "data/input.rds"
obj <- readRDS(input_file)

# Extract UMAP coordinates and cell type annotations
if (class(obj)[1] == "Seurat") {
  umap_df <- data.frame(
    UMAP_1 = obj@reductions$umap@cell.embeddings[, 1],
    UMAP_2 = obj@reductions$umap@cell.embeddings[, 2],
    celltype = obj@meta.data$X16_main
  )
} else if (class(obj)[1] == "ArchRProject") {
  umap_red <- getReducedDimensions(obj, name = "UMAP")
  umap_df <- data.frame(
    UMAP_1 = umap_red[, 1],
    UMAP_2 = umap_red[, 2],
    celltype = obj$cellColData$X16_main
  )
} else {
  stop("Input object must be Seurat or ArchRProject")
}

# Define color palette
color_id <- ...

# Create base UMAP plot
p1 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0, alpha = 1, position = "jitter") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = "white")
  ) +
  theme(
    legend.title = element_blank(),
    legend.key = element_rect(fill = 'white'),
    legend.text = element_text(size = 20, family = "sans"),
    legend.key.size = unit(1, 'cm')
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = color_id)

# Generate no-legend plot and extract legend
p2 <- p1 + NoLegend()
leg <- cowplot::get_legend(p1)
legend_plot <- cowplot::plot_grid(leg)

# Save plots (modify paths as needed)
ggsave("./Plot_of_Global_data/Global_snRNA_subtype_UMAP(Nolegend).png", p2, width = 20, height = 20, dpi = 300)
ggsave("./Plot_of_Global_data/Global_snRNA_subtype_UMAP_legend.pdf", egg::set_panel_size(legend_plot, width = unit(100, "mm"), height = unit(300, "mm")))
