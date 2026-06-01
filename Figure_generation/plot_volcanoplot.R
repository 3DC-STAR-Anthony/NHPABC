# Load required packages
library(ggplot2)
library(dplyr)
library(scales)
library(egg)

###################################################################----Cell Proportion Age Effect Volcano Plot----###################################################################
# Read input regression result table
merged_data <- read.csv("./Figure2_prop_new_level2/results/merged_results.csv", stringsAsFactors = FALSE)

# Filter target comparison stage
all_ltsr <- merged_data[merged_data$Stage == 'p_prop', ]
all_ltsr$X <- all_ltsr$Celltype
all_ltsr$Age_ltsr <- all_ltsr$LFSR
all_ltsr$Age_postmean <- all_ltsr$Estimate

# Classify regulation trend
result_df <- all_ltsr %>%
  mutate(
    trend = case_when(
      Estimate > 0.005 & Age_ltsr > 0.99 ~ "up",
      Estimate < -0.005 & Age_ltsr > 0.99 ~ "down",
      TRUE ~ "other"
    )
  )

# Create volcano plot
p_volcano <- ggplot(result_df, aes(x = Estimate, y = -log10(LFSR), color = trend)) +
  geom_point(size = 0.8, alpha = 1) +
  scale_color_manual(
    values = c("up" = "red", "down" = "blue", "other" = "gray70"),
    labels = c("up" = "Up-regulated", "down" = "Down-regulated", "other" = "Not significant")
  ) +
  # Threshold reference lines
  geom_vline(xintercept = c(-0.005, 0.005), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", alpha = 0.5) +
  labs(
    x = "Estimate (Effect Size)",
    y = "-log10(Age LFSR)",
    title = "Age Effects on Cell Types",
    subtitle = "Red: Up-regulated (Estimate > 0.005 & LFSR > 0.99)\nBlue: Down-regulated (Estimate < -0.005 & LFSR > 0.99)",
    color = "Regulation"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 6, color = "black"),
    legend.text = element_text(size = 6, color = "black"),
    legend.title = element_text(size = 6, color = "black"),
    plot.title = element_text(size = 8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  ylim(c(0, 8))

# Export figure
ggsave(
  "p_prop_neuron_LFSR_volcanoplot.pdf",
  egg::set_panel_size(p_volcano, width = unit(24, "mm"), height = unit(28, "mm")),
  width = 10, height = 10, units = "in", dpi = 300, limitsize = FALSE
)
