# Load packages
library(Seurat)
library(ggplot2)
library(egg)

###################################################################----Age Group Gene Expression Violin Plot----###################################################################
# Load multi-region Seurat objects
rds <- list()
region_list <- c("PFC","Hf","Str","Tha","Hyt","MB","PM","Cer")
for (f in region_list){
  obj_path <- paste0("./scRNA/Each_single_region/",f,".rds")
  a <- readRDS(obj_path)
  a@meta.data$Age_group <- factor(a@meta.data$Age_group, levels = c('Young','Middle age','Old','Exceptionally old'))
  rds[[f]] <- a
}

# Age group color palette
age_group_color <- c(
  'Young' = "#cc340c",
  'Middle age' = "#f18800",
  'Old' = "#9ec417",
  'Exceptionally old' = "#44c1f0"
)

# Reusable violin plot function
VlnPlot_new2 <- function(seurat_obj, subtypes, cell, gene_sel, age_group_color){
  Idents(seurat_obj) <- subtypes
  seurat_sub <- subset(seurat_obj, subtype == cell)
  
  p <- VlnPlot(
    seurat_sub,
    features = gene_sel,
    pt.size = 0,
    group.by = "Age_group",
    cols = age_group_color,
    split.by = 'Age_group',
    raster = FALSE
  ) + theme(
    axis.text.x = element_text(color = '#000000', size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.line = element_line(linewidth = 0.5),
    plot.title = element_text(size = 6, hjust = 0.5, family = "sans", color = "black", face = "italic")
  )
  return(p)
}

# Generate violin plot for NLGN1 in L5/6-NP Ex PFC
p_vln <- VlnPlot_new2(
  seurat_obj = rds[["PFC"]],
  subtypes = "subtype",
  cell = "L5/6-NP Ex",
  gene_sel = "NLGN1",
  age_group_color = age_group_color
)

# Export figure
ggsave(
  paste0("NLGN1"," ","L5-6-NP Ex PFC","_receptor","_Vlnplot.pdf"),
  egg::set_panel_size(p_vln, width = unit(16, "mm"), height = unit(19, "mm")),
  height = 5, width = 5, limitsize = FALSE
)
