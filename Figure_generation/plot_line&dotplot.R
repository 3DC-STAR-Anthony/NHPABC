# Load packages
library(ggplot2)
library(egg)

###################################################################----Gene Aging Expression Trend Line & Scatter Plot----###################################################################
# Global color palettes
region_color <- c(
  PFC = '#1B9E77', Hf = '#D95F02', Str = '#7570B3', Tha = '#E7298A',
  Hyt = '#66A61E', MB = '#E6AB02', PM = '#A6761D', Cer = '#666666'
)
label_color <- MetBrewer::metbrewer_pal("Hiroshige", n = 20)

# Reusable line plot drawing function
draw_lineplot <- function(mean_exp_sel, text_size = 6, gene = NULL, plotwidth = 15, plotheight = 25,
                          color_by = c("region", "label"), use_subtype_linetype = FALSE){
  # Select color palette
  if (color_by == "region") {
    color_palette <- region_color
    legend_name <- "Region"
  } else {
    color_palette <- label_color
    legend_name <- "Subtype"
  }
  
  # Base ggplot
  p <- ggplot(mean_exp_sel, aes(x = age, y = mean.exp, color = .data[[color_by]]))
  
  # Set line type
  if (color_by == "region" && use_subtype_linetype && "subtype" %in% colnames(mean_exp_sel)) {
    p <- p + aes(linetype = subtype)
  } else {
    p <- p + aes(linetype = NULL)
  }
  
  # Scatter point + smooth regression line
  p <- p +
    geom_point(size = 0.15) +
    geom_smooth(method = "lm", span = 0.6, se = FALSE, linewidth = 0.375, alpha = 0.8) +
    scale_color_manual(values = color_palette, name = legend_name) +
    scale_x_continuous(breaks = c(5, 10, 20, 30), labels = c("5", "10", "20", "30")) +
    theme_bw() +
    xlab("Age(years)") +
    ylab("Expression") +
    ggtitle(gene) +
    theme(
      axis.text.x = element_text(size = text_size, angle = 0, hjust = 0.5, vjust = 1, family = "sans", color = "black"),
      axis.text.y = element_text(size = text_size, family = "sans", color = "black"),
      axis.title.x = element_text(size = text_size, family = "sans", color = "black"),
      axis.title.y = element_text(size = text_size, family = "sans", color = "black"),
      legend.title = element_text(size = text_size, family = "sans", color = "black"),
      legend.text = element_text(size = text_size, family = "sans", color = "black"),
      panel.grid.major = element_line(linewidth = 0.1875, color = "grey80"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = text_size, hjust = 0.5, family = "sans", color = "black", face = "italic")
    )
  
  # Save figure
  ggsave(
    paste0(gene, "_lineplot.pdf"),
    egg::set_panel_size(p, width = unit(plotwidth, "mm"), height = unit(plotheight, "mm")),
    height = 5, width = 5, limitsize = FALSE
  )
  return(p)
}

# Load full mean expression table
mean_exp <- read.csv("./gene_age_mean_expression.csv", stringsAsFactors = FALSE)

# Filter data for EPHA4 gene in PFC region
gene_sel <- "EPHA4"
label_sel_s <- label_sel[c(1,2,3)]
mean_exp_sel <- mean_exp[mean_exp$label %in% label_sel_s & mean_exp$region == "PFC" & mean_exp$gene %in% gene_sel, ]

# Generate line & dot plot
p_epha4_trend <- draw_lineplot(
  mean_exp_sel,
  text_size = 6,
  gene = "EPHA5",
  plotwidth = 16,
  plotheight = 19,
  color_by = "label",
  use_subtype_linetype = FALSE
)
