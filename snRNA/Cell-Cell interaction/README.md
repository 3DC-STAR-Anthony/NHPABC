# Cell-cell Interaction Analysis for snRNA-seq in NHPABC
This project performs cell-cell interaction analysis on single-nuclei RNA-seq data using the [`CellChat`](https://github.com/jinworks/CellChat) package. The analysis is executed at the level of region-specific age groups, with input data formatted as `{region}_{age_group}.rds` (e.g., `PFC_Young.rds`).
# Input
- A Seurat object in `.rds` format named: `<region>.<age_group>.rds` (e.g., `PFC_Young.rds`)

- The object must contain:
  - RNA assay
  - Metadata columns including:
    - Annotion (e.g., `Subtype`)
    - Age group (e.g., `Young`)
- How to Get Assay and Meta
```r
source("~/getAllAgeGroupsCountsAndMeta.R")
result <- getAllAgeGroupsCountsAndMeta(
    obj = rds, 
    cluster = "subtype", 
    ident = "Age_group", 
    save_path = "~/01.Project/NHPABC/Figure5/01.CellChat_data/",
    prefix = "PFC"
)
```
# Creat CellChat object
```r
source("~/BuildAgegroupCellChatFromCountAndMeta.R")
cellchat <- BuildAgegroupCellChatFromCountAndMeta(
  base_path = "~/01.Project/NHPABC/Figure5/01.CellChat_data/",
  prefix = "PFC",
  save_path = "~/01.Project/NHPABC/Figure5/02.CellChat_objects/"
)
```
# CellChat Standard Process
```r
source("~/Cellchat_process_with_db_option.R")
batch_process(
    input_dir = "~/01.Project/NHPABC/Figure5/02.CellChat_objects/",
    output_dir = "~/01.Project/NHPABC/Figure5/03.Processed_CellChat/",
    prefix = "PFC",
    use_custom_db = TRUE
)
```
# Extract LR table from CellChat result
```r
source("~/extract_cellchat_communication.R")
extract_cellchat_communication(work_dir = "~/01.Project/NHPABC/Figure5/03.Processed_CellChat_prob_truncatedMean/02.result/")
```
# Get Aging-related ligand and receptor
Inter-cluster communication probabilities were computed via the computeCommunProb function (truncatedMean = 0) at each age stage. Interaction pairs with mean communication probability > 0.002 or P-value ≥ 0.05 were filtered out. An interaction was considered aging-relevant if either its ligand or receptor was included in the DEG list.

Detail process please check: 03.NHPABC_CellChat_ARLR_truncatedMean.ipynb
# Dependencies
Before running the script, make sure to install the following environment and R packages:
- R version 4.4.1 (2024-06-14)
- Seurat_5.3.0
- CellChat_2.1.2
