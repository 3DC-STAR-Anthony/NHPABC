# 将人数据中的细胞类型归成两大类神经元（兴奋性、抑制性）和四种胶质细胞（小胶质细胞，星形胶质细胞，少突胶质细胞和少突胶质前体细胞）

#Load package
library(dplyr)
library(data.table)

#######################处理Renbing的人脑数据#########################
#加载任兵的cCRE data
ren <- data.frame(fread("~/01.Project/NHPABC/Figure6/06.Conservation/cCREs.bed",nThread = 40))
head(ren,2)

#归成两大类神经元（兴奋性、抑制性）和四种胶质细胞（小胶质细胞，星形胶质细胞，少突胶质细胞和少突胶质前体细胞）
ren$celltype <- unlist(lapply(X = ren$V5, FUN = function(x) {return(strsplit(x, split = "_")[[1]][[1]])}))

ren[grepl("ASC",ren$celltype),"celltype"] <- "ASC"
ren[grepl("D12",ren$celltype),"celltype"] <- "MSN"
ren[grepl("D1",ren$celltype),"celltype"] <- "MSN"
ren[grepl("D2",ren$celltype),"celltype"] <- "MSN"
ren[grepl("COP",ren$celltype),"celltype"] <- "OPC"

message("合并后，任兵数据的celltype有")
print(unique(ren$celltype))

ren$main_celltype <- ifelse(ren$celltype %in% c("ACBGM","ASC"),"Astrocyte",
                     ifelse(ren$celltype %in% c("MGC"),"Microglia",
                     ifelse(ren$celltype %in% c("OPC"),"OPC",
                     ifelse(ren$celltype %in% c("OGC"),"Oligodendrocyte",
                     ifelse(ren$celltype %in% c("AMY","CBGRC","CT","ERC","ET","ITL23","ITL34","ITL4","ITL45","ITL5","ITL6","ITV1C","L6B","NP","PIR","SUB","TP"),"Ex",
                     ifelse(ren$celltype %in% c("BFEXA","BNGA","CBINH","CNGA","MSN","FOXP2","LAMP5","PKJ","PVALB","PV","SNCG","SST","THMGA","VIP"),"In",
                     ifelse(ren$celltype %in% c("EC","SMC"),"VS",
                                                            "Other")))))))

#拆分成各个main_celltype文件, 并利用bedtools将有重复区间的peak进行merge, 具体代码在：~/processing_conservation.sh
#保存和bedtools处理
path <- "~/01.human_maincelltype_cCRE/"
for(f in unique(ren$main_celltype)){
    a <- ren[ren$main_celltype == f, c("V1","V2","V3","main_celltype")]
    fwrite(a, file = paste0(path,f,".renbing_human_brain.bed"), row.names = F, col.names = F, sep = "\t", nThread = 40)
}

#######################处理NHPABC的数据#########################
allcCRE <- read.csv("/hwfssz3/PS_JLU/zhangxiao6/Figure6/CPM_filterr4_n4/All_59st_cCRE_list_cpm_4_individual_n_4.csv")





