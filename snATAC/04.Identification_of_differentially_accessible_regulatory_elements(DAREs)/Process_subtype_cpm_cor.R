# ============================================================================
# Correlation Analysis with Permutation Test
# ============================================================================
cor.perm <- function(x, y, nperm = 10000) {
  r.obs <- cor(x = x, y = y)
  P.par <- cor.test(x = x, y = y)$p.value
  r.per <- sapply(1:nperm, function(i) cor(x = x, y = sample(y)))
  r.per <- c(r.per, r.obs)
  P.per <- sum(abs(r.per) >= abs(r.obs)) / (nperm + 1)
  return(list(r.obs = r.obs, P.par = P.par, P.per = P.per))
}

# ============================================================================
# Process Subtype Correlation from RDS File
# ============================================================================
Process_subtype_cpm_cor <- function(input_rds, output_dir, mc_cores = 140) {
  library(parallel)
  
  # Sample information
  monkey_age <- data.frame(
    monkey = c('Y1','Y2','Y3','Y4','Y5','Y6','Y7','M1','M2','M3','M4','M5',
               'O1','O2','O3','O4','O5','O6','V1','V2','V3','V4','V5','V6'),
    age = c(5,6,5,5,5,5,5,10,11,12,10,10,22,23,22,22,22,23,29,29,28,31,29,31)
  )
  rownames(monkey_age) <- monkey_age$monkey
  
  # Process file
  subtype_name <- tools::file_path_sans_ext(basename(input_rds))
  message("Processing: ", subtype_name, " ------ ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  cpm_matrix <- readRDS(input_rds)
  message("  Matrix: ", nrow(cpm_matrix), " features x ", ncol(cpm_matrix), " samples")
  
  if (is.null(colnames(cpm_matrix))) stop("CPM matrix must have column names!")
  
  # Match samples
  common_samples <- intersect(colnames(cpm_matrix), rownames(monkey_age))
  if (length(common_samples) < ncol(cpm_matrix)) {
    warning("Some samples not found, using ", length(common_samples), " matched samples")
    cpm_matrix <- cpm_matrix[, common_samples, drop = FALSE]
  }
  
  # Calculate correlations
  message("  Calculating correlations for ", nrow(cpm_matrix), " features...")
  results <- mclapply(seq_len(nrow(cpm_matrix)), function(i) {
    cor_result <- cor.perm(x = as.numeric(cpm_matrix[i,]), y = monkey_age[colnames(cpm_matrix), "age"])
    data.frame(
      peaks = rownames(cpm_matrix)[i],
      cor = cor_result$r.obs,
      pvalue = cor_result$P.par,
      perm = cor_result$P.per
    )
  }, mc.cores = mc_cores) %>% do.call(rbind, .)
  
  # Save full results
  full_path <- file.path(output_dir, paste0(subtype_name, "_cor_results.rds"))
  saveRDS(results, file = full_path)
  message("  Full results saved: ", full_path)
  
  # Filter and save
  filtered <- results[abs(results$cor) > 0.5 & results$pvalue < 0.05 & results$perm < 0.1, ]
  filtered_path <- file.path(output_dir, paste0(subtype_name, "_cor_filtered.rds"))
  saveRDS(filtered, file = filtered_path)
  message("  Filtered results saved: ", filtered_path, " (", nrow(filtered), " features, ", 
          round(nrow(filtered)/nrow(results)*100, 2), "%)")
  
  message("Completed: ", subtype_name, " ------ ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  return(filtered)
}
