# **Analysis of cell type abundance across aging in NHPABC**
To ensure robust detection, we retained cell types only if they are represented across at least two age groups and ≥2 individuals (each with ≥10 cells)
## Calculate proportions
```r
source("~/Calculate_proportions.R")
meta <- read.csv("~/PFC_All_neuorn_meta_data.csv")
df <- calculate_proportions(meta)
```
## Identify progressive linear alterations
```r
source("~/Calculate_LTSR.R")
analysis_results <- analyze_age_effects(df)
write.csv(analysis_results$results, "PFC_All_neuron_age_effects_with_lfsr_p_prop_add_mixed.csv"), row.names = FALSE)
```
## Identify stage-specific alterations
```r
unique(df$Age_group)
copairson_group <- list(
                       'young_vs_middle'  = c('Young','Middle age'),
                       'middle_vs_old'  = c('Middle age', 'Old'),
                       'old_vs_expold'  = c('Old', 'Exceptionally old'))

for (f in seq(1:3)){
df <- df[df$Age_group %in% copairson_group[[f]],]
df$Age[df$Age_group == copairson_group[[f]][1]] <- 1
df$Age[df$Age_group == copairson_group[[f]][2]] <- 2
analysis_results <- analyze_age_effects(df)

results_1 <-  analysis_results$results
results_1$Region <- unique(df$Region)
results_1$Stage <- names(copairson_group)[f]

write.csv(results_1, paste0(outdir, unique(df$Region),"_age_effects_with_ltsr_",names(copairson_group)[f],"_s_prop_add_mixed_neuron.csv"), row.names = FALSE)
    }

```
