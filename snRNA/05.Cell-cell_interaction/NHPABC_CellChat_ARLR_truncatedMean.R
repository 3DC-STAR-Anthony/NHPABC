# ============================================================================
# NHPABC CellChat Analysis: Aging-Related Ligand-Receptor (ARLR) Identification
# Figure: Fig5 (Cell-cell communication network)
# Type: network
# Description:
#   This script processes CellChat cell-cell interaction (CCI) results across
#   eight brain regions (PFC, Hf, Str, Tha, Hyt, MB, PM, Cer) and four age
#   groups (Young, Middle_age, Old, Exceptionally_old). It merges interaction
#   data, filters low-quality events, and identifies aging-related ligand-
#   receptor pairs (ARLR) by integrating differentially expressed genes (DEGs)
#   from each aging stage.
#
#   Input: CellChat probability matrices per region and age group
#   Output: Filtered CCI table with ARLR annotations
# ============================================================================

# ---- 1. Load packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(data.table)
  library(parallel)
  library(matrixStats)
})

message("Packages loaded successfully.")

# ---- 2. Set paths and parameters ----
# Input/output directories (modify as needed)
work_dir        <- "path/to/your/working/directory/"
input_dir       <- file.path(work_dir, "CellChat_prob_results")
output_dir      <- file.path(work_dir, "output")
custom_lr_db    <- file.path(work_dir, "interaction_input_CellChatDB_custom.csv")
complex_ref_file<- file.path(work_dir, "CellChatDB.human.complex.gene_list.rds")
deg_file        <- file.path(work_dir, "DEG_all_filter.csv")

# Brain regions and age groups
regions    <- c("PFC", "Hf", "Str", "Tha", "Hyt", "MB", "PM", "Cer")
age_groups <- c("Young", "Middle_age", "Old", "Exceptionally_old")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- 3. Load reference data ----
# Custom ligand-receptor interaction database
interaction_df <- read.csv(custom_lr_path, row.names = 1)

# Complex gene reference (for handling multi-subunit receptors/ligands)
complex_ref <- readRDS(complex_ref_file)

message("Reference data loaded. LR database: ", nrow(interaction_df), " interactions.")

# ---- 4. Define merge function ----
#' Merge CellChat results across brain regions and age groups
#'
#' @param input_dir Directory containing CellChat probability CSV files
#' @param output_dir Directory to save merged results
#' @param regions Vector of brain region names
#' @param age_groups Vector of age group names
#' @param interaction_df Custom LR interaction database
#' @param complex_ref Complex gene reference list
#' @param n_threads Number of threads for fread/fwrite
#'
#' @return List containing per-region results and combined dataset
merge_cellchat_precise_v2 <- function(input_dir, output_dir,
                                      regions = c("PFC", "Hf", "Str", "Tha", "Hyt", "MB", "PM", "Cer"),
                                      age_groups = c("Young", "Middle_age", "Old", "Exceptionally_old"),
                                      interaction_df = interaction_df,
                                      complex_ref = complex_ref,
                                      n_threads = 8) {

  # Internal function: process a single region
  .process_region <- function(region) {
    message("\nProcessing region: ", region, " ------ ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

    # 1. Read all age group files for this region
    age_data <- list()
    for (age in age_groups) {
      file_pattern <- paste0(region, "_", age, "_subsetCommunication_prob_truncatedMean.csv")
      file_path <- file.path(input_dir, file_pattern)

      if (file.exists(file_path)) {
        df <- fread(file_path, data.table = FALSE, nThread = n_threads)

        # Map interaction annotations from custom database
        df$annotation <- interaction_df[match(df$interaction_name, interaction_df$interaction_name), "annotation"]
        df <- df[!is.na(df$annotation), ]

        df$interaction_name_2 <- interaction_df[match(df$interaction_name_2, interaction_df$interaction_name_2), "annotation"]
        df <- df[!is.na(df$annotation), ]

        df$pathway_name <- interaction_df[match(df$interaction_name, interaction_df$interaction_name), "pathway_name"]

        # Resolve complex ligand/receptor genes
        inter_df <- df
        tmp <- merge(inter_df, complex_ref, by.x = 3, by.y = 2, all.x = TRUE)
        inter_df <- merge(tmp, complex_ref, by.x = 4, by.y = 2, all.x = TRUE)
        inter_df$gene_name.x[is.na(inter_df$gene_name.x)] <- inter_df$ligand[is.na(inter_df$gene_name.x)]
        inter_df$gene_name.y[is.na(inter_df$gene_name.y)] <- inter_df$receptor[is.na(inter_df$gene_name.y)]
        colnames(inter_df)[c(12, 13)] <- c("ligand_gene", "receptor_gene")
        df <- inter_df

        # Add metadata and create unique interaction ID
        df$region <- region
        df$age_group <- gsub("_", " ", age)
        df$interaction_id <- paste0(
          df$source, "#", df$target, "#",
          df$ligand, "#", df$receptor, "#",
          df$interaction_name, "#", df$interaction_name_2, "#",
          df$pathway_name, "#", df$annotation, "#",
          df$region, "#", df$ligand_gene, "#", df$receptor_gene
        )

        # Remove duplicates
        df %>% distinct(interaction_id, .keep_all = TRUE) -> dfun
        rownames(dfun) <- dfun$interaction_id

        age_data[[age]] <- dfun
        message("  ", nrow(dfun), " unique interactions from ", age)
      } else {
        warning("File not found: ", file_path)
        age_data[[age]] <- NULL
      }
    }

    # 2. Merge all age groups and find unique interactions
    merged_data <- Reduce(rbind, age_data)
    unique_interactions <- unique(merged_data$interaction_id)
    message("\n  Total unique interactions across age groups: ", length(unique_interactions))

    # 3. Create complete template with all unique interactions
    merged_data %>% distinct(interaction_id, .keep_all = TRUE) -> template
    if (nrow(template) == length(unique_interactions)) {
      template$prob <- 0
      template$pval <- NA
      rownames(template) <- template$interaction_id
    } else {
      stop("unique_interactions number does not match template length!")
    }

    # 4. Fill in data for each age group
    complete_age_data <- list()
    for (age in names(age_data)) {
      df_old <- age_data[[age]]
      df_new <- template
      df_new$age_group <- age
      df_new[rownames(df_old), "prob"] <- df_old$prob
      df_new[rownames(df_old), "pval"] <- df_old$pval
      complete_age_data[[age]] <- df_new
    }

    # 5. Combine all age groups for this region
    region_complete <- Reduce(rbind, complete_age_data)
    message("  Final interaction count: ", nrow(region_complete))

    return(region_complete)
  }

  # Process all regions
  final_results <- list()
  for (region in regions) {
    region_result <- .process_region(region)
    if (!is.null(region_result)) {
      final_results[[region]] <- region_result

      # Save per-region results
      filename <- file.path(output_dir, paste0(region, "_merge_cellchat_truncatedMean_LR_pair_result.csv"))
      fwrite(region_result, filename, nThread = n_threads)
      message("  Saved: ", basename(filename))
    }
  }

  # Combine all regions
  message("\nCombining results from all regions...")
  all_regions_complete <- Reduce(rbind, final_results)
  message("Total interactions across all regions: ", nrow(all_regions_complete))

  # Save combined results
  saveRDS(all_regions_complete, file.path(output_dir, "01.NHPABC_CellChat_truncatedMean_merge.rds"))
  fwrite(all_regions_complete, file.path(output_dir, "01.NHPABC_CellChat_truncatedMean_merge.csv"), nThread = n_threads)

  return(list(
    Each_region_results = final_results,
    all_regions_complete = all_regions_complete
  ))
}

# ---- 5. Execute data merging ----
# NOTE: This step is time-consuming (estimated: 20-60 mins depending on hardware).
# If already run, load the merged result directly (see Step 6).

merge_all <- merge_cellchat_precise_v2(
  input_dir  = input_dir,
  output_dir = output_dir,
  regions    = regions,
  age_groups = age_groups,
  interaction_df = interaction_df,
  complex_ref    = complex_ref,
  n_threads  = 8
)

# Alternatively, load pre-computed merged results:
# merge_all <- fread(file.path(output_dir, "01.NHPABC_CellChat_truncatedMean_merge.csv"),
#                    nThread = 8, data.table = FALSE)

# ---- 6. Prepare merged data ----
cci <- merge_all  # or merge_all$all_regions_complete if from RDS

# ---- 7. Cell type grouping (Source -> Target) ----
# Classify source and target cells into broad categories: Neuron, Glia, Other

cci <- cci %>%
  mutate(source_celltype = case_when(
    grepl("Ast|Mic|ODC|OPC", source) ~ "Glia",
    grepl("VS|Ependymal|Fibroblast|BEC|Choroid|Macrophage|Pericyte|SMC|UBC", source) ~ "Other",
    TRUE ~ "Neuron"
  ))

cci <- cci %>%
  mutate(target_celltype = case_when(
    grepl("Ast|Mic|ODC|OPC", target) ~ "Glia",
    grepl("VS|Ependymal|Fibroblast|BEC|Choroid|Macrophage|Pericyte|SMC|UBC", target) ~ "Other",
    TRUE ~ "Neuron"
  ))

cci$CellToCell <- paste0(cci$source_celltype, "_", cci$target_celltype)

# ---- 8. Integrate data across age groups ----
# Core function: standardize CCI data across age groups for cross-stage comparison

inter_all <- cci
inter_all$ID <- paste0(inter_all$interaction_id, "#", inter_all$CellToCell)
unique_ID <- distinct(inter_all, ID, .keep_all = TRUE)$ID

message("Total CCI entries: ", nrow(inter_all))
message("Unique interaction IDs: ", length(unique_ID), " (expected: ~1/4 of total entries)")

# Split by age group and extract probability/pvalue columns
df.list <- list()
for (ag in unique(inter_all$age_group)) {
  message("Processing age group: ", ag)
  df <- subset(inter_all, age_group == ag)[, c(3, 4, 1, 2, 12, 13, 7:10, 14, 15, 16, 19, 20, 5:6)]
  prob_col <- gsub("_", " ", ag)
  colnames(df)[16:17] <- c(prob_col, paste0("pval_", gsub(" ", "_", ag)))
  rownames(df) <- df$ID
  df <- df[unique_ID, ]
  df.list[[ag]] <- df
}

# Merge age groups: Exceptionally_old | Middle_age | Old | Young
# Order: Young -> Middle_age -> Old -> Exceptionally_old
age_order <- c("Young", "Middle_age", "Old", "Exceptionally_old")
df_ordered <- df.list[age_order]

inter_single <- cbind(
  df_ordered[[4]],           # Exceptionally_old (base columns)
  df_ordered[[2]][, 16:17],  # Middle_age
  df_ordered[[3]][, 16:17],  # Old
  df_ordered[[1]][, 16:17]   # Young
)

# Reorder columns: metadata + probs + pvals
inter_single <- inter_single[, c(1:15, 16, 18, 20, 22, 17, 19, 21, 23)]
rownames(inter_single) <- NULL

# ---- 9. Initial filtering on CCI data ----
# Multi-step filtering to retain high-quality, biologically meaningful events:
#   1. Remove P-value > 0.05 (all four age groups)
#   2. Remove entries with missing P-values
#   3. Remove low-probability events (mean prob across ages < 0.002)

# Helper: print filtering statistics
print_filter_stats <- function(df_before, df_after, step_name = "Filtering") {
  n_before <- nrow(df_before)
  n_after  <- nrow(df_after)
  n_filtered <- n_before - n_after
  filter_ratio <- ifelse(n_before == 0, 0, n_filtered / n_before * 100)

  message(sprintf("%s before: %d", step_name, n_before))
  message(sprintf("%s after:  %d", step_name, n_after))
  message(sprintf("%s removed: %d (%.2f%%)", step_name, n_filtered, filter_ratio))
}

# Step 1: Filter P-value > 0.05
inter_single1 <- inter_single[which(
  rowMaxs(as.matrix(inter_single[, c(20:23)]), na.rm = TRUE) < 0.05
), ]
print_filter_stats(inter_single, inter_single1, step_name = "Step1_Pvalue")

# Step 2: Filter missing P-values
inter_single1$missing_value <- ifelse(
  is.na(inter_single1$pval_Young) |
  is.na(inter_single1$pval_Middle_age) |
  is.na(inter_single1$pval_Old) |
  is.na(inter_single1$pval_Exceptionally_old),
  "TRUE", "FALSE"
)
inter_single2 <- subset(inter_single1, missing_value == "FALSE")
print_filter_stats(inter_single1, inter_single2, step_name = "Step2_Missing")

# Step 3: Filter low-probability events
mean_probs <- rowMeans2(as.matrix(inter_single2[, 16:19]))
inter_single3 <- inter_single2[which(mean_probs > 0.002), ]
print_filter_stats(inter_single2, inter_single3, step_name = "Step3_LowProb")

# Summary of probability distribution
sorted_means <- sort(mean_probs, decreasing = TRUE)
message("\nProbability distribution (before Step 3):")
message("  95th percentile: ", sorted_means[as.integer(length(sorted_means) * 0.95)])
message("  50th percentile: ", sorted_means[as.integer(length(sorted_means) * 0.50)])

# ---- 10. Identify aging-related ligand-receptor pairs (ARLR) ----
# Core: annotate LR pairs where ligand or receptor is a DEG at any aging stage

# Load DEG results
deg <- read.csv(deg_file)
deg <- deg[, -c(1, 2, 3)]  # Remove index columns
colnames(deg)[4:8] <- c("fdr", "celltype", "trend", "region", "stage")
deg$label <- paste0(deg$subtype, "_", deg$region)
deg <- deg[deg$subtype != "VS", ]

# Split DEGs by aging stage
pdeg     <- deg[deg$stage == "pDEG", ]
Early    <- deg[deg$stage == "Early", ]
Late     <- deg[deg$stage == "Late", ]
Very_Late <- deg[deg$stage == "Very_Late", ]

# Create source/target labels for matching with DEGs
inter_single3$source_label <- paste0(inter_single3$source, "_", inter_single3$region)
inter_single3$target_label <- paste0(inter_single3$target, "_", inter_single3$region)

# Remove cell types not present in DEG data
deg_label <- unique(deg$label)
cci_label <- unique(c(inter_single3$source_label, inter_single3$target_label))
outcelltype <- setdiff(cci_label, deg_label)

inter_single3$TF <- ifelse(
  inter_single3$source_label %in% outcelltype |
  inter_single3$target_label %in% outcelltype, "T", "F"
)
inter_single4 <- subset(inter_single3, TF == "F")[, -27]

# ---- 11. Annotate aging-related ligands and receptors ----
#' Annotate LR pairs with aging-related genes from DEG data
#'
#' @param cci_data Filtered CCI data frame
#' @param deg_data_list Named list of DEG data frames per stage
#'
#' @return Named list of annotated CCI data per stage
annotate_aging_lr <- function(cci_data, deg_data_list) {
  cci_list <- list()

  for (stage in names(deg_data_list)) {
    cci_tmp <- cci_data
    cci_tmp$ligand_age   <- NA_character_
    cci_tmp$receptor_age <- NA_character_
    deg_stage <- deg_data_list[[stage]]

    # Annotate ligand: if ligand_gene is a DEG in source cell type
    for (sl in unique(cci_tmp$source_label)) {
      genes <- deg_stage$gene[deg_stage$label == sl]
      idx <- cci_tmp$source_label == sl
      cci_tmp$ligand_age[idx] <- ifelse(
        cci_tmp$ligand_gene[idx] %in% genes,
        cci_tmp$ligand_gene[idx], NA
      )
    }

    # Annotate receptor: if receptor_gene is a DEG in target cell type
    for (tl in unique(cci_tmp$target_label)) {
      genes <- deg_stage$gene[deg_stage$label == tl]
      idx <- cci_tmp$target_label == tl
      cci_tmp$receptor_age[idx] <- ifelse(
        cci_tmp$receptor_gene[idx] %in% genes,
        cci_tmp$receptor_gene[idx], NA
      )
    }

    cci_list[[stage]] <- cci_tmp
  }

  return(cci_list)
}

# Run annotation for all aging stages
deg_list <- list(
  pDEG      = pdeg,
  Early     = Early,
  Late      = Late,
  Very_Late = Very_Late
)

CCI_list <- annotate_aging_lr(inter_single4, deg_list)

# ---- 12. Compile ARLR results ----
# Add stage labels and merge
for (stage in names(CCI_list)) {
  CCI_list[[stage]]$stage <- stage
}

merge_all <- Reduce(rbind, CCI_list)

# Define ARLR: TRUE if either ligand or receptor is aging-related
merge_all$ARLR <- !is.na(merge_all$ligand_age) | !is.na(merge_all$receptor_age)

# Extract ARLR-positive entries
merge_all_ARLR <- merge_all[merge_all$ARLR, ]

# Remove duplicates (same ID may match multiple DEGs)
merge_all_un <- merge_all_ARLR %>% distinct(ID, .keep_all = TRUE)

message("\nARLR summary:")
message("  Total annotated entries: ", nrow(merge_all))
message("  ARLR-positive (with duplicates): ", nrow(merge_all_ARLR))
message("  ARLR-unique (deduplicated): ", nrow(merge_all_un))

# ---- 13. Save results ----
fwrite(merge_all,     file.path(output_dir, "03.NHPABC_CellChat_truncatedMean_final_all.csv"),     nThread = 8)
fwrite(merge_all_ARLR, file.path(output_dir, "04.NHPABC_CellChat_truncatedMean_final_ARLR.csv"),  nThread = 8)
fwrite(merge_all_un,   file.path(output_dir, "04.NHPABC_CellChat_truncatedMean_final_ARLR_unique.csv"), nThread = 8)

message("\nAll results saved to: ", output_dir)
message("Done! ------ ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
