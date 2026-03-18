# Load CellChat library
#library(CellChat)

# Built-in path for custom database files
custom_lr_path <- "~/interaction_input_CellChatDB_update_20250909_LGY_modfiy_V1.csv"
custom_gene_path <- "~/geneInfo_input_CellChatDB_update_20250911_LGY_modfiy_V1.csv"

# Core function to process CellChat object
process_cellchat <- function(input_cellchat_path, output_path, use_custom_db = TRUE) {
    # Load CellChat object from RDS file
    cellchat <- readRDS(input_cellchat_path)
    
    # Initialize default CellChat database
    CellChatDB <- CellChatDB.human
    
    # Select database based on parameter
    if (use_custom_db) {
        # Load custom ligand-receptor and gene information database
        # Built-in path for custom database files
        options(stringsAsFactors = FALSE)
        CellChatDB$interaction <- read.csv(custom_lr_path, row.names = 1)
        CellChatDB$geneInfo <- read.csv(custom_gene_path, row.names = 1)
        message("Using custom database")
    } else {
        message("Using original CellChatDB database")
    }
    cellchat@DB <- CellChatDB
    
    # Run CellChat analysis pipeline
    source("~/identifyOverExpressedGenes_chenge.R")
    cellchat <- subsetData(cellchat) # Subset the expression data of signaling genes to save computation cost
    cellchat <- identifyOverExpressedGenes_chenge(cellchat,thresh.p = 1)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE,type = "truncatedMean",trim = 0)#10% trimming 
    cellchat <- filterCommunication(cellchat, min.cells = 10)#Filter out cell groups with fewer than 100 cells
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    # Save processed CellChat object
    saveRDS(cellchat, output_path)
    message("Saved to: ", output_path)
    return(cellchat)
}

# Batch process all age groups
batch_process <- function(input_dir, output_dir, prefix, use_custom_db = TRUE) {
    age_groups <- c("Old", "Exceptionally old", "Middle age", "Young")
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Process each age group sequentially
    lapply(age_groups, function(group) {
        # Sanitize group name (replace spaces with underscores) for file naming
        sanitized <- gsub(" ", "_", group)
        input_path <- file.path(input_dir, paste0(prefix, "_", sanitized, "_cellchat.rds"))
        # Add database type suffix to output filename
        db_suffix <- ifelse(use_custom_db, "customDB", "originalDB")
        output_path <- file.path(output_dir, paste0(prefix, "_", sanitized, "_", db_suffix, "_processed.rds"))
        
        # Process current age group
        process_cellchat(input_path, output_path, use_custom_db)
    })
}

# Usage examples
# 1. Use custom database (default)
#batch_process(
#    input_dir = "~/01.Project/NHPABC/Figure5/02.CellChat_objects/",
#    output_dir = "~/01.Project/NHPABC/Figure5/03.Processed_CellChat/",
#    prefix = "PFC",
#    use_custom_db = TRUE
#)

# 2. Use original CellChatDB database
# batch_process(
#     input_dir = "~/01.Project/NHPABC/Figure5/02.CellChat_objects/",
#     output_dir = "~/01.Project/NHPABC/Figure5/03.Processed_CellChat/",
#     prefix = "PFC",
#     use_custom_db = FALSE
# )
