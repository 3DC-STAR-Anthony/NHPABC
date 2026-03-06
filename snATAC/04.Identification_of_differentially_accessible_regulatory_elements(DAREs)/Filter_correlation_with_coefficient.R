# ============================================================================
# Filter Correlation Results with Coefficient Match
# ============================================================================
# Description: This function performs secondary filtering by matching correlation
#              results with coefficient analysis results. It compares correlation
#              trends with coefficient trends and filters features that match
#              in direction and have significant coefficient values.

filter_correlation_with_coefficient <- function(subtype_name, cor_results, coeff_results, output_dir, coef_threshold = 0.005) {
  message("Filtering ", subtype_name, " ------ ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  # Merge coefficient information
  cor_results$coef <- coeff_results[match(cor_results$peaks, coeff_results$gene), "age.logFC"]
  cor_results$coef_p <- coeff_results[match(cor_results$peaks, coeff_results$gene), "age.H_p"]
  
  # Determine trends
  cor_results$cor_trend <- ifelse(cor_results$cor > 0, "Up", "Down")
  cor_results$coef_trend <- ifelse(cor_results$coef > 0, "Up", "Down")
  cor_results$match <- ifelse(cor_results$cor_trend == cor_results$coef_trend, "match", "Notmatch")
  
  # Save raw merged results
  saveRDS(cor_results, file.path(output_dir, paste0(subtype_name, "_cor_addcoeff_raw.rds")))
  
  # Filter by trend match and coefficient threshold
  matched_results <- cor_results[cor_results$match == "match" & abs(cor_results$coef) > coef_threshold, ]
  
  # Save filtered results
  saveRDS(matched_results, file.path(output_dir, paste0(subtype_name, "_cor_addcoeff_filter.rds")))
  
  message(subtype_name, " - pDARE: ", nrow(matched_results))
  return(matched_results)
}
