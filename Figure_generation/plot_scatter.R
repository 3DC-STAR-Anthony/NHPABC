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

###################################################################----Fig1C----###################################################################

# ===================== Load data =====================
# Load Seurat scRNA RDS, replace path with your rds path
scrna_obj <- readRDS("./scrna_seurat.rds")
meta_rna <- scrna_obj@meta.data %>%
  filter(Tech == "RNA")

# If you need ATAC data plot, uncomment below and replace ArchR path
# library(ArchR)
# atac_archr <- loadArchRProject("./ArchRProject/")
# meta_atac <- atac_archr@cellColData %>%
#   as.data.frame() %>%
#   filter(Tech == "ATAC")

# ===================== Plot setting =====================
# Region color mapping
Region_color <- c(
  PFC = '#1B9E77',
  Hf = '#D95F02',
  Str = '#7570B3',
  Tha = '#E7298A',
  Hyt = '#66A61E',
  MB = '#E6AB02',
  PM = '#A6761D',
  Cer = '#666666'
)
# Cell type & Region factor levels
Cell_level <- c(
  'Cortical Ex',
  'CGE/MGE-derived Inh',
  'Hippocampal CA Ex',
  'Hippocampal DG Ex',
  'SPN',
  'MB-derived Inh',
  'TPH_N',
  'Thalamic Ex',
  'Mixed neuron',
  'Astrocyte',
  'ODC',
  'Microglia',
  'OPC',
  'Ependymal',
  'Choroid plexus',
  'VS',
  'MLI',
  'Cerebellar GC'
)
Region_level <- c('PFC','Hf','Str','Tha','Hyt','MB','PM','Cer')
numeric_labels <- seq_along(Cell_level)

# ===================== RNA data plot =====================
df_rna <- as.data.frame(table(meta_rna$cell_type, meta_rna$Region), stringsAsFactors = F)
df_stat_rna <- df_rna %>%
  group_by(Var1) %>%
  summarise(Region = Var2, Prop = Freq/sum(Freq)*100) %>%
  mutate(Var1 = factor(Var1, Cell_level),
         Region = factor(Region, Region_level))

p_rna <- ggplot(df_stat_rna, aes(x=Var1, y=Region)) +
  geom_point(aes(size = Prop, color = Region)) +
  theme_cowplot() +
  scale_size_continuous(range = c(0.5,2)) +
  scale_color_manual(values = Region_color) +
  scale_x_discrete(labels = numeric_labels) +
  scale_y_discrete() +
  ylab("Region") +
  xlab("") +
  theme(
    text = element_text(size = 6, family = "sans", color = "black"),
    axis.text.x = element_text(size = 6, family = "sans", color = "black", vjust = 0.5, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family="sans",size=6, colour="black",hjust = 0.5,vjust = 0.5),
    axis.text.y=element_text(family="sans",size=6, colour="black"),
    axis.line = element_line(linewidth = 0.25),
    axis.line.y = element_line(linewidth = 0.25),
    axis.ticks.x = element_line(linewidth = 0.25),
    axis.ticks.y = element_line(linewidth = 0.25),
    legend.position = "none"
  )

# Save RNA plot
ggsave('Global_scRNA_subtype_distribution_no_legend.pdf',
       egg::set_panel_size(p_rna, width=unit(80, "mm"), height=unit(20, "mm")),
       height = 10, width = 10, limitsize = FALSE)

# ===================== ATAC data plot (uncomment to run) =====================
# df_atac <- as.data.frame(table(meta_atac$cell_type, meta_atac$Region), stringsAsFactors = F)
# df_stat_atac <- df_atac %>%
#   group_by(Var1) %>%
#   summarise(Region = Var2, Prop = Freq/sum(Freq)*100) %>%
#   mutate(Var1 = factor(Var1, Cell_level),
#          Region = factor(Region, Region_level))
#
# p_atac <- ggplot(df_stat_atac, aes(x=Var1, y=Region)) +
#   geom_point(aes(size = Prop, color = Region)) +
#   theme_cowplot() +
#   scale_size_continuous(range = c(0.5,2)) +
#   scale_color_manual(values = Region_color) +
#   scale_x_discrete(labels = numeric_labels) +
#   scale_y_discrete() +
#   ylab("Region") +
#   xlab("") +
#   theme(
#     text = element_text(size = 6, family = "sans", color = "black"),
#     axis.text.x = element_text(size = 6, family = "sans", color = "black", vjust = 0.5, hjust = 0.5),
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(family="sans",size=6, colour="black",hjust = 0.5,vjust = 0.5),
#     axis.text.y=element_text(family="sans",size=6, colour="black"),
#     axis.line = element_line(linewidth = 0.25),
#     axis.line.y = element_line(linewidth = 0.25),
#     axis.ticks.x = element_line(linewidth = 0.25),
#     axis.ticks.y = element_line(linewidth = 0.25),
#     legend.position = "none"
#   )
#
# ggsave('Global_snATAC_subtype_distribution_no_legend.pdf',
#        egg::set_panel_size(p_atac, width=unit(80, "mm"), height=unit(20, "mm")),
#        height = 10, width = 10, limitsize = FALSE)

###################################################################----Fig1C----###################################################################








