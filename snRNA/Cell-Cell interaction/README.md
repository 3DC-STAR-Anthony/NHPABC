# Cell-cell Interaction Analysis for snRNA-seq in NHPABC
This project performs cell-cell interaction analysis on single-nuclei RNA-seq data using the [`CellChat`](https://github.com/jinworks/CellChat) package. The analysis is executed at the level of region-specific age groups, with input data formatted as `{region}_{age_group}.rds` (e.g., `PFC_Young.rds`).
# Input
- A Seurat object in `.rds` format named: `<region>.<age_group>.rds` (e.g., `PFC_Young.rds`)

- The object must contain:
  - RNA assay
  - Metadata columns including:
    - Annotion (`Subtype`)
    - Age_group (`Young`)
- How to Get Assay and Meta
```r
result <- getAllAgeGroupsCountsAndMeta(
    obj = rds, 
    cluster = "subtype", 
    ident = "Age_group", 
    save_path = "~/01.Project/NHPABC/Figure5/01.CellChat_data/",
    prefix = "PFC"
)
