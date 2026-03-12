library(ArchR)
library(dplyr)
library(tidyr)

analyze_age_qc_correlation <- function(proj, metric = "TSS", output_path, fdr_method = "fdr", fdr_threshold = 0.05) {
  if (!metric %in% c("TSS", "FRIP")) {
    stop("metric must be either 'TSS' or 'FRIP'")
  }
  
  metric_col <- ifelse(metric == "TSS", "TSSEnrichment", "FRIP")
  meta <- as.data.frame(proj@cellColData)
  
  if (!metric_col %in% colnames(meta)) {
    stop("Column '", metric_col, "' not found in cell metadata")
  }
  if (!"subtype" %in% colnames(meta)) {
    stop("Column 'subtype' not found in cell metadata")
  }
  if (!"Age" %in% colnames(meta)) {
    stop("Column 'Age' not found in cell metadata")
  }
  
  age_lm <- meta %>%
    group_by(subtype) %>%
    summarise(
      lm_model = list(lm(as.formula(paste(metric_col, "~ Age")), na.action = na.omit)),
      Intercept = coef(lm_model[[1]])[1],
      Age_Coefficient = coef(lm_model[[1]])[2],
      Age_Pvalue = summary(lm_model[[1]])$coefficients[2, 4],
      R_squared = summary(lm_model[[1]])$r.squared,
      n_cells = n()
    ) %>%
    arrange(desc(abs(Age_Coefficient)))
  
  age_lm <- as.data.frame(age_lm)
  age_lm$Age_FDR <- p.adjust(age_lm$Age_Pvalue, method = fdr_method)
  
  if (any(age_lm$Age_FDR == 0)) {
    non_zero_fdr <- age_lm$Age_FDR[age_lm$Age_FDR != 0]
    if (length(non_zero_fdr) > 0) {
      min_fdr <- min(non_zero_fdr)
      age_lm$Age_FDR[age_lm$Age_FDR == 0] <- min_fdr
    }
  }
  
  age_lm$neg_log10_FDR <- -log10(age_lm$Age_FDR)
  neg_log10_threshold <- -log10(fdr_threshold)
  age_lm$significant <- age_lm$Age_FDR < fdr_threshold
  
  age_summary <- age_lm %>%
    select(-lm_model) %>%
    as.data.frame()
  
  write.csv(age_summary, output_path, row.names = FALSE)
  
  return(age_summary)
}
