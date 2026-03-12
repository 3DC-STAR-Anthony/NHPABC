library(dplyr)

# Calculate proportions and statistics
    calculate_proportions <- function(df) {
        result <- df %>%
            group_by(Sample, Age, Age_group, Subtype) %>%
            summarise(n = n(), .groups = 'drop') %>%
            group_by(Sample, Age) %>%
            mutate(Proportion = n / sum(n),
                   Total_cell = sum(n),
                   CPM = round(Proportion, 6) * 1e6) %>%
            ungroup()
        
        cat("Dimensions after proportion calculation:", dim(result), "\n")
        
        # Mark samples with n >= 10 instead of filtering them out
        result <- result %>%
            mutate(is_valid_sample = n >= 10)
        
        # Filter criteria: 
        # 1. Each cell type must have at least 2 valid samples (n >= 10) in each age group
        # 2. Each cell type must be present in at least 4 age groups
        valid_Subtypes <- result %>%
            filter(is_valid_sample) %>%
            group_by(Subtype, Age_group) %>%
            summarise(valid_sample_count = n_distinct(Sample), .groups = 'drop') %>%
            group_by(Subtype) %>%
            filter(valid_sample_count >= 2) %>%
            summarise(valid_age_groups = n_distinct(Age_group), .groups = 'drop') %>%
            filter(valid_age_groups >= 4) %>%
            pull(Subtype)
        
        cat("Number of eligible cell types:", length(valid_Subtypes), "\n")
        cat("Eligible cell types:", paste(valid_Subtypes, collapse = ", "), "\n")
        
        # Keep all data for eligible cell types (including samples with n < 10)
        result <- result %>%
            filter(Subtype %in% valid_Subtypes)
        
        cat("Dimensions after filtering:", dim(result), "\n")
        
        return(result)
    }
