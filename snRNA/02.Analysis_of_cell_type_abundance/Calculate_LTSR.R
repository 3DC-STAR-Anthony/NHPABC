# Load required packages
library(lme4)
library(lmerTest)
library(MASS)
library(multcomp)  # For multiple comparison correction
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Data preprocessing and checking
data_preprocessing <- function(meta) {
  # Ensure correct data types
  meta$Sample <- as.factor(meta$Sample)
  meta$Subtype <- as.factor(meta$Subtype)
  
  # Check for zero or negative Total_cell values
  if(any(meta$Total_cell <= 0)) {
    warning("Some Total_cell values are <= 0, removing these rows")
    meta <- meta[meta$Total_cell > 0, ]
  }
  
  # Check sample counts per cell type
  Subtype_counts <- table(meta$Subtype)
  print("Sample counts per Subtype:")
  print(Subtype_counts)
  
  # Filter out cell types with too few samples (less than 3)
  valid_Subtypes <- names(Subtype_counts[Subtype_counts >= 3])
  meta <- meta[meta$Subtype %in% valid_Subtypes, ]
  
  return(meta)
}

# 2. Model fitting function - with convergence checking
fit_Subtype_model <- function(ct_data, Subtype_name) {
  # Check if data is sufficient
  if(nrow(ct_data) < 3) {
    cat("Insufficient data, skipping cell type:", Subtype_name, "\n")
    return(NULL)
  }
  
  # Check data variability
  if(sd(ct_data$n) == 0) {
    cat("No variability in cell counts, skipping cell type:", Subtype_name, "\n")
    return(NULL)
  }
  
  # Use only glmer model
  model <- NULL
  tryCatch({
    # Try Poisson mixed model
    model <- glmer(n ~ Age + Batch + (1 | Sample) + offset(log(Total_cell)),
                   family = poisson(link = "log"), 
                   data = ct_data,
                   control = glmerControl(optimizer = "bobyqa", 
                                        optCtrl = list(maxfun = 2e5)))
    
    # Check if model converged
    if(!is.null(model)) {
      # Check convergence status
      conv_info <- model@optinfo$conv
      if(!is.null(conv_info$opt) && conv_info$opt != 0) {
        cat("Model did not converge, skipping cell type:", Subtype_name, "\n")
        cat("Convergence info:", conv_info$opt, "-", conv_info$messages, "\n")
        return(NULL)
      }
    }
    
  }, error = function(e) {
    cat("Model fitting failed, skipping cell type:", Subtype_name, "-", e$message, "\n")
    return(NULL)
  })
  
  return(model)
}

# 3. ltsr calculation function
calculate_ltsr <- function(model) {
  if (is.null(model)) {
    return(NA)
  }
  
  tryCatch({
    # Get fixed effect coefficients and standard errors
    summ <- summary(model)
    coef_table <- summ$coefficients
    
    # Calculate z-values
    z_values <- coef_table[, "Estimate"] / coef_table[, "Std. Error"]
    
    # Calculate two-tailed p-values
    p_values <- 2 * pnorm(-abs(z_values))
    
    # Calculate ltsr (Local False Sign Rate)
    # ltsr = 1 - posterior probability (correct sign)
    # Use two-tailed p-values to calculate ltsr
    ltsr <- pmin(p_values, 1 - p_values)
    
    return(ltsr)
  }, error = function(e) {
    cat("ltsr calculation failed:", e$message, "\n")
    return(NA)
  })
}

# 4. Main analysis function - with empty result handling
analyze_age_effects <- function(meta) {
  # Data preprocessing
  meta_clean <- data_preprocessing(meta)
  Subtypes <- unique(meta_clean$Subtype)
  
  # Store results
  results <- data.frame()
  models_list <- list()
  ltsr_list <- list()
  
  # Analyze each cell type
  for(ct in Subtypes) {
    cat("Analyzing cell type:", ct, "\n")
    
    ct_data <- meta_clean[meta_clean$Subtype == ct, ]
    
    # Fit model
    model <- fit_Subtype_model(ct_data, ct)
    
    if(!is.null(model)) {
      # Extract coefficients and statistics
      tryCatch({
        summ <- summary(model)
        coef_table <- summ$coefficients
        
        # Calculate ltsr
        ltsr_values <- calculate_ltsr(model)
        ltsr_list[[ct]] <- ltsr_values
        
        # Extract Age effect results
        if("Age" %in% rownames(coef_table)) {
          age_coef <- coef_table["Age", ]
          
          result_row <- data.frame(
            Subtype = ct,
            Estimate = age_coef["Estimate"],
            StdError = age_coef["Std. Error"],
            Zvalue = age_coef["z value"],
            Pvalue = age_coef["Pr(>|z|)"],
            ltsr = ifelse("Age" %in% names(ltsr_values), ltsr_values["Age"], NA),
            N_samples = nrow(ct_data),
            Model_type = "glmer",
            stringsAsFactors = FALSE
          )
          
          results <- rbind(results, result_row)
          models_list[[ct]] <- model
        }
      }, error = function(e) {
        cat("Coefficient extraction failed, skipping cell type:", ct, "-", e$message, "\n")
      })
    }
  }
  
  # Multiple testing correction (only if results are not empty)
  if(nrow(results) > 0) {
    results$Padj <- p.adjust(results$Pvalue, method = "fdr")
    results$Significant <- results$Padj < 0.05
    results$ltsr_significant <- ifelse(!is.na(results$ltsr), results$ltsr < 0.05, FALSE)
  } else {
    cat("Warning: No cell types were successfully analyzed\n")
  }
  
  return(list(results = results, models = models_list, ltsr = ltsr_list))
}
