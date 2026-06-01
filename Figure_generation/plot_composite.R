# Load package
library(ggplot2)

# Read cCRE peak list
cpm_list <- read.csv("./All_cCRE_list.csv", stringsAsFactors = FALSE)

# Clean peak chr format
cpm_list$peaks <- gsub("-", ":", cpm_list$peaks)
cpm_list$peaks <- gsub("NC_", "NC-", cpm_list$peaks)

# Filter valid gene-cCRE pairs
genePeakPair.corr_filter <- subset(
  genePeakPair.df,
  TSS_distance <= c(500000) &
    genePeakName %in% unique(cpm_list$peaks) &
    estimate > c(0.3)
)

# Generate random shuffled background pairs
genePeakPair.rdmShuf <- subset(genePeakPair.df, class == "rdmShuf")
genePeakPair.rdmShuf$pair <- paste0(genePeakPair.rdmShuf$gene, "_", genePeakPair.rdmShuf$name)
genePeakPair.rdmShuf <- genePeakPair.rdmShuf[sample(nrow(genePeakPair.rdmShuf), c(3329941)), ]

# Subset valid correlation & shuffle data
genePeakPair.corr_f <- genePeakPair.corr_filter[, c(16,17,18,19)]
genePeakPair.rdmShuf_f <- genePeakPair.rdmShuf[, c(16,17,18,19)]

# Merge real correlation and random shuffle data
genePeakPair.corr_shuf <- rbind(genePeakPair.rdmShuf_f, genePeakPair.corr_f)

# Plot composite density distribution
p_density <- ggplot(genePeakPair.corr_shuf, aes(x = estimate)) +
  geom_density(aes(fill = class)) +
  scale_fill_manual(values = c("corr" = "#FF3333", "rdmShuf" = "grey")) +
  geom_vline(xintercept = c(0.3), color = "red", linetype = "solid", linewidth = c(0.5), alpha = c(1)) +
  theme_bw()

# Export figure
ggsave(
  "gene_cCRE_correlation_density.pdf",
  p_density,
  width = c(10), height = c(6), limitsize = FALSE
)
