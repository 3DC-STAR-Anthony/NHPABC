#library(CellChat)

# 内置自定义数据库路径
custom_lr_path <- "/hwfssz3/PS_JLU/laiguangyao/02.Database/01.Construct_LR_database/04.CellChat_v2/interaction_input_CellChatDB_update_20250909_LGY_modfiy_V1.csv"
custom_gene_path <- "/hwfssz3/PS_JLU/laiguangyao/02.Database/01.Construct_LR_database/04.CellChat_v2/geneInfo_input_CellChatDB_update_20250911_LGY_modfiy_V1.csv"

# 处理CellChat对象的核心函数
process_cellchat <- function(input_cellchat_path, output_path, use_custom_db = TRUE) {
    # 加载CellChat对象
    cellchat <- readRDS(input_cellchat_path)
    
    # 初始化数据库
    CellChatDB <- CellChatDB.human
    
    # 根据参数选择数据库
    if (use_custom_db) {
        # 加载自定义数据库
        # 内置自定义数据库路径
        options(stringsAsFactors = FALSE)
        CellChatDB$interaction <- read.csv(custom_lr_path, row.names = 1)
        CellChatDB$geneInfo <- read.csv(custom_gene_path, row.names = 1)
        message("使用自定义数据库")
    } else {
        message("使用原始CellChatDB数据库")
    }
    cellchat@DB <- CellChatDB
    
    # 运行分析流程
    cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
    #cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedGenes_chenge(cellchat,thresh.p = 1)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    #cellchat <- projectData(cellchat, PPI.human)
    #cellchat <- computeCommunProb(cellchat, raw.use = TRUE,type = "triMean",trim = 0.1)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE,type = "truncatedMean",trim = 0)#10% ½Ø¶ÏµÄÆ½¾ùÖµ
    cellchat <- filterCommunication(cellchat, min.cells = 0)#filter out cell number less then 100
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    # 保存结果
    saveRDS(cellchat, output_path)
    message("已保存至: ", output_path)
    return(cellchat)
}

# 批量处理所有年龄组
batch_process <- function(input_dir, output_dir, prefix, use_custom_db = TRUE) {
    age_groups <- c("Old", "Exceptionally old", "Middle age", "Young")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    lapply(age_groups, function(group) {
        sanitized <- gsub(" ", "_", group)
        input_path <- file.path(input_dir, paste0(prefix, "_", sanitized, "_cellchat.rds"))
        db_suffix <- ifelse(use_custom_db, "customDB", "originalDB")
        output_path <- file.path(output_dir, paste0(prefix, "_", sanitized, "_", db_suffix, "_processed.rds"))
        
        process_cellchat(input_path, output_path, use_custom_db)
    })
}

# 使用示例
# 1. 使用自定义数据库（默认）
#batch_process(
#    input_dir = "/hwfssz3/PS_JLU/laiguangyao/01.Project/NHPABC/Figure5/02.CellChat_objects/",
#    output_dir = "/hwfssz3/PS_JLU/laiguangyao/01.Project/NHPABC/Figure5/03.Processed_CellChat/",
#    prefix = "PFC",
#    use_custom_db = TRUE
#)

# 2. 使用原始数据库
# batch_process(
#     input_dir = "/hwfssz3/PS_JLU/laiguangyao/01.Project/NHPABC/Figure5/02.CellChat_objects/",
#     output_dir = "/hwfssz3/PS_JLU/laiguangyao/01.Project/NHPABC/Figure5/03.Processed_CellChat/",
#     prefix = "PFC",
#     use_custom_db = FALSE
# )

CellChatDB <- CellChatDB.human
#interaction_input <- read.csv("/mnt/11/monkey_aging_brain/figure4/09.Aging_related_CellChat/01.Built_LR_database/Merge_CellChat_neuronChat_interactionDB.csv",row.names = 1)
#interaction_input <- read.csv("/hwfssz3/PS_JLU/laiguangyao/02.Database/01.Construct_LR_database/01.LR_pair_database/ALL_Tech_LRDB_Deduplication_final.csv",row.names = 1)
#geneIfo <- read.csv("/hwfssz3/PS_JLU/laiguangyao/02.Database/01.Construct_LR_database/03.Other_Version/geneInfo_input_CellChatDB_update.csv",row.names = 1)
#CellChatDB$interaction <- interaction_input
#CellChatDB$geneInfo <- geneIfo
#cellchat@DB <- CellChatDB
#cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
#cellchat <- identifyOverExpressedGenes(cellchat)
#source("/hwfssz1/CS_CELL/cs_cell/laiguangyao/script/01.snRNA/Cell-cell_interaction/CellChat/identifyOverExpressedGenes_chenge.R")
#cellchat <- identifyOverExpressedGenes_chenge(cellchat)
#cellchat <- identifyOverExpressedGenes_chenge(cellchat,thresh.p = 1)
#cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
#cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "triMean", trim = 0.25) ##
#cellchat <- computeCommunProb(cellchat, raw.use = TRUE,type = "truncatedMean",trim = 0)#10% ½Ø¶ÏµÄÆ½¾ùÖµ
#cellchat <- filterCommunication(cellchat, min.cells = 10) #filter out cell number less then 100
#cellchat <- computeCommunProbPathway(cellchat)
    