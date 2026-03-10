# Load packages
library(dplyr)
library(data.table)

####################### Process Ren's human brain data #########################
# Load Ren's cCRE data
ren <- data.frame(fread("~/01.Project/NHPABC/Figure6/06.Conservation/cCREs.bed", nThread = 40))
head(ren, 2)

# Assign to two major neuron categories and four glial cell types
ren$celltype <- unlist(lapply(X = ren$V5, FUN = function(x) {return(strsplit(x, split = "_")[[1]][[1]])}))

# Standardize cell type names
ren[grepl("ASC", ren$celltype), "celltype"] <- "ASC"
ren[grepl("D12", ren$celltype), "celltype"] <- "MSN"
ren[grepl("D1", ren$celltype), "celltype"] <- "MSN"
ren[grepl("D2", ren$celltype), "celltype"] <- "MSN"
ren[grepl("COP", ren$celltype), "celltype"] <- "OPC"

message("After merging, cell types in Ren's data:")
print(unique(ren$celltype))

# Assign to main cell type categories
ren$main_celltype <- ifelse(ren$celltype %in% c("ACBGM", "ASC"), "Astrocyte",
                     ifelse(ren$celltype %in% c("MGC"), "Microglia",
                     ifelse(ren$celltype %in% c("OPC"), "OPC",
                     ifelse(ren$celltype %in% c("OGC"), "Oligodendrocyte",
                     ifelse(ren$celltype %in% c("AMY", "CBGRC", "CT", "ERC", "ET", "ITL23", "ITL34", "ITL4", "ITL45", "ITL5", "ITL6", "ITV1C", "L6B", "NP", "PIR", "SUB", "TP"), "Ex",
                     ifelse(ren$celltype %in% c("BFEXA", "BNGA", "CBINH", "CNGA", "MSN", "FOXP2", "LAMP5", "PKJ", "PVALB", "PV", "SNCG", "SST", "THMGA", "VIP"), "In",
                     ifelse(ren$celltype %in% c("EC", "SMC"), "VS",
                                                            "Other")))))))

# Split into separate files for each main_celltype, then use bedtools to merge overlapping peaks
# (specific code in: ~/processing_conservation.sh)
# Save and process with bedtools
path <- "~/01.human_maincelltype_cCRE/"
for(f in unique(ren$main_celltype)){
    a <- ren[ren$main_celltype == f, c("V1", "V2", "V3", "main_celltype")]
    fwrite(a, file = paste0(path, f, ".renbing_human_brain.bed"), row.names = F, col.names = F, sep = "\t", nThread = 40)
}

####################### Process NHPABC data #########################
allcCRE <- read.csv("~/All_subtype_cCRE.csv")
allcCRE <- allcCRE %>%
  mutate(bigcluster = case_when(
      grepl("Ast", subtype) ~ "Ast",
      grepl("Mic|Macrophage", subtype) ~ "Mic",
      grepl("ODC", subtype) ~ "ODC",
      grepl("OPC|COP", subtype) ~ "OPC",
      grepl(" Ex", subtype) ~ "Ex",
      grepl("rhombic lip", subtype) ~ "Ex",
      grepl("Cerebellar", subtype) ~ "Ex",
      grepl(" Inh", subtype) ~ "Inh",
      grepl("SPN", subtype) ~ "Inh",
      grepl("MLI|PLI|Purkinje", subtype) ~ "Inh",
      grepl("A10 DAn|CHO|Splatter neuron|A9 DAn", subtype) ~ "Other neuron",
      grepl("Ependymal", subtype) ~ "Ependymal",
      grepl("VS", subtype) ~ "VS",
      TRUE ~ "Unknown"
  ))

# Filter to keep only relevant cell types
allcCRE_sub <- allcCRE[allcCRE$bigcluster %in% c('Ast', 'Ex', 'Mic', 'Inh', 'ODC', 'OPC', 'VS'), ]

# Save each subtype to separate bed files
setwd("~/01.Project/NHPABC/Figure6/06.Conservation/02.NHPABC_cCRE_bed")
for(f in unique(allcCRE_sub$subtype)){
    a <- allcCRE_sub[allcCRE_sub$subtype == f, c("seqnames", "start", "end", "peaks", "bigcluster")]
    bigcluster <- unique(a$bigcluster)
    fwrite(a[, 1:4], file = paste0(f, ".", bigcluster, ".nhpabc.bed"), row.names = F, col.names = F, sep = "\t", nThread = 40)
}
