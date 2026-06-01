# Load required packages
library(ggplot2)
library(dplyr)
library(scales)
library(egg)

###################################################################----DARE/DEG Up/Down Frequency Distribution Histogram----###################################################################
# Read full dataset
pdare_all <- read.csv("./pDARE_all_result.csv", stringsAsFactors = FALSE)

# ---------------- Down-regulated pDARE Histogram ----------------
# Count peak frequency table for Down trend
pdare_st_down <- as.data.frame(table(pdare_all[pdare_all$cor_trend == "Down", ]$peaks), stringsAsFactors = FALSE)
pdare_st_down <- subset(pdare_st_down, Freq != 0)

# Calculate statistic ratio
n1 <- length(pdare_st_down[pdare_st_down$Freq == "1", ]$Var1)
sum47 <- sum(pdare_st_down[pdare_st_down$Freq >= 47, ]$Freq)

# Draw histogram for Down pDARE
p_down_hist <- ggplot(pdare_st_down, aes(x = Freq)) +
  geom_bar(stat = "count", position = "stack", fill = "#619399") +
  geom_vline(xintercept = 46.5, color = "black", linetype = "dashed", linewidth = 0.5) +
  xlab("") +
  ylab("No. of pDAREs") +
  ggtitle("") +
  scale_x_continuous(limits = c(0, 94), breaks = c(1, 30, 60, 90), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::comma_format(), expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 1, family = "sans", color = "black"),
    axis.text.y = element_text(size = 6, family = "sans", color = "black"),
    axis.title.x = element_text(size = 6, family = "sans", color = "black"),
    axis.title.y = element_text(size = 6, family = "sans", color = "black"),
    legend.title = element_text(size = 6, family = "sans", color = "black"),
    legend.text = element_text(size = 6, family = "sans", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 8, hjust = 0.5, family = "sans", color = "black")
  )

# Export Down histogram
ggsave(
  "NHPABC_pDAREs_Down_Distribution_hist.pdf",
  egg::set_panel_size(p_down_hist, width = unit(35, "mm"), height = unit(14.3, "mm")),
  height = 4, width = 15, limitsize = FALSE
)

# ---------------- Up-regulated pDARE Histogram ----------------
# Count peak frequency table for Up trend
pdare_st_up <- as.data.frame(table(pdare_all[pdare_all$cor_trend == "Up", ]$peaks), stringsAsFactors = FALSE)
pdare_st_up <- subset(pdare_st_up, Freq != 0)

# Calculate statistic ratio
n1_up <- length(pdare_st_up[pdare_st_up$Freq == "1", ]$Var1)
sum47_up <- sum(pdare_st_up[pdare_st_up$Freq >= 47, ]$Freq)

# Draw histogram for Up pDARE
p_up_hist <- ggplot(pdare_st_up, aes(x = Freq)) +
  geom_bar(stat = "count", position = "stack", fill = "#FEBB44") +
  geom_vline(xintercept = 46.5, color = "black", linetype = "dashed", linewidth = 0.5) +
  xlab("") +
  ylab("No. of pDAREs") +
  ggtitle("") +
  scale_x_continuous(limits = c(0, 94), breaks = c(1, 30, 60, 90), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::comma_format(), expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 1, family = "sans", color = "black"),
    axis.text.y = element_text(size = 6, family = "sans", color = "black"),
    axis.title.x = element_text(size = 6, family = "sans", color = "black"),
    axis.title.y = element_text(size = 6, family = "sans", color = "black"),
    legend.title = element_text(size = 6, family = "sans", color = "black"),
    legend.text = element_text(size = 6, family = "sans", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 8, hjust = 0.5, family = "sans", color = "black")
  )

# Export Up histogram
ggsave(
  "NHPABC_pDAREs_Up_Distribution_hist.pdf",
  egg::set_panel_size(p_up_hist, width = unit(35, "mm"), height = unit(14.3, "mm")),
  height = 4, width = 15, limitsize = FALSE
)
