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

# Load data
scrna_obj <- readRDS("./scrna_seurat.rds")
meta_rna <- scrna_obj@meta.data %>% filter(Tech == "RNA")

# Plot settings
Region_color <- c(
  PFC = '#1B9E77', Hf = '#D95F02', Str = '#7570B3', Tha = '#E7298A',
  Hyt = '#66A61E', MB = '#E6AB02', PM = '#A6761D', Cer = '#666666'
)
Cell_level <- c('Cortical Ex','CGE/MGE-derived Inh','Hippocampal CA Ex','Hippocampal DG Ex','SPN','MB-derived Inh','TPH1+','Thalamic Ex','Splatter neuron','Astrocyte','ODC','Microglia','OPC','Ependymal','Choroid plexus','VS','MLI','Cerebellar GC')
Region_level <- c('PFC','Hf','Str','Tha','Hyt','MB','PM','Cer')
numeric_labels <- seq_along(Cell_level)

# RNA plot
df_rna <- as.data.frame(table(meta_rna$cell_type, meta_rna$Region), stringsAsFactors=FALSE)
df_stat_rna <- df_rna %>% group_by(Var1) %>% summarise(Region=Var2, Prop=Freq/sum(Freq)*100) %>% mutate(Var1=factor(Var1, Cell_level), Region=factor(Region, Region_level))

p_rna <- ggplot(df_stat_rna, aes(x=Var1, y=Region)) +
  geom_point(aes(size=Prop, color=Region)) +
  theme_cowplot() +
  scale_size_continuous(range=c(0.5,2)) +
  scale_color_manual(values=Region_color) +
  scale_x_discrete(labels=numeric_labels) +
  ylab("Region") + xlab("") +
  theme(
    text=element_text(size=6, family="sans", color="black"),
    axis.text.x=element_text(size=6, family="sans", color="black", vjust=0.5, hjust=0.5),
    axis.title.y=element_text(family="sans", size=6, color="black", hjust=0.5, vjust=0.5),
    axis.text.y=element_text(family="sans", size=6, color="black"),
    axis.line=element_line(linewidth=0.25), axis.line.y=element_line(linewidth=0.25),
    axis.ticks.x=element_line(linewidth=0.25), axis.ticks.y=element_line(linewidth=0.25),
    legend.position="none"
  )
ggsave('XXXX_subtype_distribution_no_legend.pdf', egg::set_panel_size(p_rna, width=unit(80,"mm"), height=unit(20,"mm")), height=10, width=10, limitsize=FALSE)

###################################################################----Other Plot----###################################################################

# Load data
df_early <- read.csv("./early_count.csv", stringsAsFactors=FALSE)
df_late <- read.csv("./late_count.csv", stringsAsFactors=FALSE)
df_verylate <- read.csv("./verylate_count.csv", stringsAsFactors=FALSE)
df_early$stage <- "Early"
df_late$stage <- "Late"
df_verylate$stage <- "Very late"
sdare_df <- bind_rows(df_early, df_late, df_verylate)

# Plot settings
stage_col <- c(Early="#FF6A6A", Late="#4169E1", "Very late"="#008B8B")

# Plot
plot1 <- ggplot(sdare_df, aes(x=Var1, y=Freq, group=stage, color=stage)) +
  geom_point(size=0.2) +
  coord_flip() +
  scale_color_manual(values=stage_col) +
  scale_y_reverse(labels=scales::comma_format(), expand=c(0,0), limits=c(20000,0)) +
  scale_x_discrete(position="top") +
  labs(x="", y="No. of pDAREs") +
  theme_bw() +
  theme(
    axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5, size=6, family="sans", color="black"),
    axis.title.y=element_text(size=6, family="sans", color="black"),
    axis.text.y=element_text(size=6, family="sans", color="black"),
    legend.title=element_text(size=6, color="black"),
    legend.text=element_text(size=6, color="black"),
    legend.position="left",
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  )
ggsave("XXXX_number_dotplot.pdf", egg::set_panel_size(plot1, width=unit(166,"mm"), height=unit(15,"mm")), height=4, width=15, limitsize=FALSE)
