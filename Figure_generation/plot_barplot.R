# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(future)
library(showtext)
showtext_auto()
library(cowplot)
library(egg)
library(MetBrewer)
library(ggrepel)

###################################################################----Cell Type Proportion Stacked Barplot----###################################################################
# Load Seurat object
scrna_obj <- readRDS("./scrna_seurat.rds")
meta <- scrna_obj@meta.data

# Global multi-group color palettes, switch palette according to fill variable

fill_col <- XXXX #Please add this information!

# Factor order preset
cluster_lev <- levels(meta$cluster2)

# Calculate cell count matrix
df_region_total <- as.data.frame(table(meta$Region), stringsAsFactors = FALSE)
df_raw <- as.data.frame(table(meta$Region, meta$cluster2), stringsAsFactors = FALSE)

# Calculate total cells per region and cell proportion
df_raw$total <- df_region_total$Freq[match(df_raw$Var1, df_region_total$Var1)]
df_raw$Per <- df_raw$Freq / df_raw$total

# Set display order of factors
df_raw$Var1 <- factor(df_raw$Var1, rev(region_lev))
df_raw$Var2 <- factor(df_raw$Var2, rev(cluster_lev))

# ------------------- Adjustable fill & color segment -------------------
# Switch fill variable between Var1(Region) / Var2(celltype), match corresponding palette
fill_var <- "Var1"
fill_pal <- region_col
# fill_var <- "Var2"
# fill_pal <- celltype_col
# ------------------------------------------------------------------------

# Build stacked percentage barplot
p_stack_bar <- ggplot(df_raw, aes(x = Var2, y = Freq)) +
  geom_bar(aes(fill = .data[[fill_var]]), stat = "identity", position = "fill") +
  scale_fill_manual(values = fill_pal) +
  theme_bw() +
  theme(
    text = element_text(size = 6, color = "black"),
    axis.text = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 6, color = "black"),
    legend.text = element_text(size = 6, color = "black"),
    legend.title = element_text(size = 6, color = "black"),
    axis.text.x = element_text(family = "sans", size = 6, angle = 45, hjust = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_flip()

# Save output PDF
ggsave(
  "XXX_barplot_distribution.pdf",
  egg::set_panel_size(p_stack_bar, width = unit(22, "mm"), height = unit(38, "mm")),
  height = 10, width = 16, limitsize = FALSE
)

###################################################################----2. pDEG/pDARE Count Stacked Barplot (DEG/DARE Count Type)----###################################################################
# process DEG/DARE table
df <- as.data.frame(table(df$label, pdare_all$coef_trend), stringsAsFactors = FALSE)
# Set cell type factor level
df$Var1 <- factor(df$Var1, lev$label)
# Trend color palette
trend_color <- c(
  Up = "#FEBB44",
  Down = "#619399"
)

# Plot stacked bar for up/down regulated pDARE counts
p_dare_bar <- ggplot(pdare_df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(
    labels = scales::comma_format(),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = trend_color, name = "Trend") +
  theme_bw() +
  labs(x = "", y = "No. of pDAREs/pDEGs", title = "") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3, family = "sans", color = "Black"),
    axis.title.y = element_text(size = 6, family = "sans", color = "Black"),
    axis.text.y = element_text(size = 6, family = "sans", color = "Black"),
    legend.title = element_text(size = 6, color = "Black"),
    legend.text = element_text(size = 6, color = "Black"),
    legend.position = "left",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Export PDF
ggsave(
  "NHPABC_pDAREs/pDEGs_subtype_number_barplot.pdf",
  egg::set_panel_size(p_dare_bar, width = unit(166, "mm"), height = unit(15, "mm")),
  height = 4, width = 15, limitsize = FALSE
)

###################################################################----GO Term Horizontal Barplot----###################################################################
# Load GO enrichment excel table
go <- readxl::read_excel("./go_enrich_result.xlsx")

# Format GO description text: capitalize first letter
go$Description <- ifelse(
  nchar(go$Description) > 0,
  paste0(
    toupper(substr(go$Description, 1, 1)),
    substr(go$Description, 2, nchar(go$Description))
  ),
  go$Description
)

# Custom color palette for GO terms, fill your color mapping below
Description_col <- c() #Please add this information!

# Set GO term display order
go$Description <- factor(go$Description, rev(unique(go$Description)))

# Plot GO enrichment barplot
p_go_bar <- ggplot(go, aes(x = -LogP_Gene, y = Description, fill = Description)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = Description_col, name = "Functional enrichment") +
  theme_classic() +
  xlab("-Log10(Q)") +
  ylab("") +
  ggtitle("Functional enrichment of ARLR DEGs") +
  theme(
    axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 1, family = "sans", color = "black"),
    axis.text.y = element_text(size = 6, family = "sans", color = "black"),
    axis.title.x = element_text(size = 6, family = "sans", color = "black"),
    legend.title = element_text(size = 6, family = "sans", color = "black"),
    legend.text = element_text(size = 6, family = "sans", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 8, hjust = 0.5, family = "sans", color = "black")
  )

# Export PDF figure
ggsave(
  "NHPABC_GOterm_barplot.pdf",
  egg::set_panel_size(p_go_bar, width = unit(12, "mm"), height = unit(46, "mm")),
  height = 10, width = 10, limitsize = FALSE
)

###################################################################----Genome Liftover Percentage Single Barplot----###################################################################
# Load background and cCRE bed data
bg_row <- fread("./Mac_GCF_037993035_raw.bed")
bg_hg38 <- fread("./Mac_GCF_03799303_acToHg38.bed")
ccre_row <- fread("./NHPABC_cCRE_raw.bed")
ccre_hg38 <- fread("./NHPABC_cCRE_MacToHg38.bed")

# Calculate liftover percentage
df_global <- data.frame(
  type = c("Background", "cCREs in macaque"),
  percentage = c(
    nrow(bg_hg38)/nrow(bg_row)*100,
    nrow(ccre_hg38)/nrow(ccre_row)*100
  )
)

# Plot parameters & color palette
lift_color <- c(
  "Background" = "#4D4D4D",
  "cCREs in macaque" = "#00BFC4"
)

# Generate bar plot
p_lift_bar <- ggplot(df_global, aes(x = type, y = percentage, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = lift_color, name = "") +
  ylab("Percentage") +
  xlab("") +
  ggtitle("Liftover to human genome") +
  ylim(0,100) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, family = "sans", color = "black"),
    axis.text.y = element_text(size = 6, family = "sans", color = "black"),
    axis.title.x = element_text(size = 6, family = "sans", color = "black"),
    axis.title.y = element_text(size = 6, family = "sans", color = "black"),
    legend.title = element_text(size = 6, family = "sans", color = "black"),
    legend.text = element_text(size = 6, family = "sans", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 8, hjust = 0.5, family = "sans", color = "black"),
    strip.text = element_text(size = 6)
  )

# Export output PDF
ggsave(
  paste0("NHPABC_Liftover_to_human_genome.pdf"),
  egg::set_panel_size(p_lift_bar, width = unit(10, "mm"), height = unit(26, "mm")),
  height = 10, width = 10, limitsize = FALSE
)
