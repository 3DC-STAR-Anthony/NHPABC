# Load ArchR package
library(ArchR)

###################################################################----ArchR Motif Footprinting Full Pipeline & Plot----###################################################################
# Set parallel threads
addArchRThreads(threads = 100)

# Load main ArchR Project
proj <- readRDS("./ArchR/output/total_xyc/total_ATAC_peakmatrix.rds")

# Load motif PWM database
pwm_2022 <- readRDS("./JASPAR2022_PWM_motif.rds")

# Link motif gene name with ArchR PWM name
jasper2024_list <- read.csv("./jasper2024_motif_gene_list.csv", stringsAsFactors = FALSE)
annotation <- data.frame(
  feature = unique(jasper2024_list$feature),
  ArchR_name = names(pwm_2022)
)
write.csv(annotation, "./jasper2022_motif_ArchR_list.csv", row.names = FALSE)

# Run motif annotation in ArchR project
proj@projectMetadata$outputDirectory <- "./ArchR/output/total_xyc/"
proj <- addMotifAnnotations(
  ArchRProj = proj,
  name = "Motif",
  motifPWMs = pwm_2022,
  force = FALSE
)

# Save updated ArchR project with motif data
saveRDS(proj, "./Figure6/10.Footprinting/total_ATAC_peakmatrix_addmotif.rds")

# ---------------------- ArchR native Footprint Plot ----------------------
# Example 1: Plot single motif footprint
p_foot_single <- plotFootprints(
  ArchRProj = proj,
  motifName = "FOXF2_1",
  groupBy = "label2",
  plotName = "FOXF2_Footprint",
  outputDirectory = "./Figure6/10.Footprinting/Footprint_Plots/",
  pal = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"),
  flank = c(250, 250),
  smoothWindow = 20,
  baseSize = 6
)
