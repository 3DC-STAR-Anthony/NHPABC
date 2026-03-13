# Figure7 MBA PFC deg recalculate

library(Seurat)
library(MAST)
library(dplyr)
library(tidyr)
library(Matrix)
library(RUVSeq)
library(edgeR)
library(parallel)
library(dplyr)





ALL_rds <- readRDS('/hwfssz3/PS_JLU/zhangxiao6/Figure7/MBA_PFC_filter_gene_reclustered.rds')
DefaultAssay(ALL_rds) <- 'RNA'
ALL_rds$Age_group[ALL_rds$Age_group == 'Middle age'] = 'Middle'
ALL_rds$Age_group[ALL_rds$Age_group == 'Exceptionally old'] = 'Very_old'
ALL_rds$Age[ALL_rds$Age == '10'] = 11
ALL_rds$Age[ALL_rds$Age == '11'] = 11
ALL_rds$Age[ALL_rds$Age == '12'] = 11
ALL_rds$Age[ALL_rds$Age == '23'] = 22
ALL_rds$Age[ALL_rds$Age == '22'] = 22
ALL_rds$Age[ALL_rds$Age== '29'] = 29
ALL_rds$Age[ALL_rds$Age== '28'] = 29
ALL_rds$Age[ALL_rds$Age== '31'] = 29
ALL_rds$cluster1 <- ALL_rds$subtype
celltypes <- unique(ALL_rds$cluster1)
ALL_rds$label <- paste(ALL_rds$cluster1, ALL_rds$Region, sep = '_')
deg_list <- read.csv('/hwfssz3/PS_JLU/zhangxiao6/Figure3_DEG/Revise_DEG_LIST/Revised_DEG_all_filter.csv')
label_list <- intersect(unique(deg_list$label), unique(ALL_rds$label))
label_list
length(label_list)
method = 'RUV_pDEG'

make.tform = function(x, norm=FALSE, u=NULL){
    if (is.null(u)){u = unique(x)}
    tmat = sapply(u, function(y){1 * (x == y)})
    if (norm){
        tmat = sweep(tmat, 2, apply(tmat, 2, sum),  '/')
    }
    return(tmat) 
}

subsetMatrixForDE = function(mat, pathdf, pctcells, pctcut=NULL, selgenes=NULL){
    if (is.null(selgenes)){
        if (!is.null(pctcut)){
            keep.genes = names(pctcells)[pctcells > pctcut]
        }
    } else {
        keep.genes = selgenes[selgenes %in% rownames(mat)]
    }
    mat = mat[keep.genes, pathdf$barcode]
    print(paste("[STATUS] Subsetting matrix to", 
                paste0(dim(mat), collapse = ' x '),'(g x c)'))
    return(mat)
}


asform = function(x){ as.formula(paste0(x, collapse='')) }

est.endtime <- function(start, last, done, total){
    curr = proc.time() # Get current time
    # Calculate how much time has passed:
    last.step = (proc.time() - last)[3]
    elapsed = (proc.time() - start)[3]
    cat(paste0(done,'/', total, '\t'))
    cat(paste0(round(last.step,1),'s\t'))
    cat(paste0(round(elapsed,1),'s\t'))
    # Estimate the remaining time:
    # Estimate the end time:
    each.time = elapsed / done
    est.left = (total - done) * each.time
    cat(paste0('Left: ', round(est.left,1),'s\t'))
    fin.time = Sys.time() + est.left
    cat(paste0('Est. ', format(fin.time, "%H:%M:%S"),'\n'))
    return(curr)
}

subsetMatrixForDE = function(mat, pathdf, pctcells, pctcut=NULL, selgenes=NULL){
    if (is.null(selgenes)){
        if (!is.null(pctcut)){
            keep.genes = names(pctcells)[pctcells > pctcut]
        }
    } else {
        keep.genes = selgenes[selgenes %in% rownames(mat)]
    }
    mat = mat[keep.genes, pathdf$barcode]
    print(paste("[STATUS] Subsetting matrix to", 
                paste0(dim(mat), collapse = ' x '),'(g x c)'))
    return(mat)
}

runRUVpsbulk = function(pathdf, mat, path,nruv){ 
    pids = as.character(unique(pathdf$Sample))
    tform = make.tform(pathdf$Sample, u=pids)
    data_ind = mat %*% tform 
    indvar = 'Sample'
    uqcols = c(indvar,path)
    dform = asform(c('~',path))
    uqobs = unique(pathdf[,uqcols])
    rownames(uqobs) = uqobs[[indvar]]
    uqobs = uqobs[colnames(data_ind),]
    design = model.matrix(dform, data=uqobs)
    d_e = DGEList(data_ind, genes = rownames(data_ind))
    keep = rowSums(cpm(d_e) > 1) >= 3
    d_e = d_e[keep, , keep.lib.sizes = FALSE]
    d_e = calcNormFactors(d_e, method = "TMM")
    d_e = estimateGLMCommonDisp(d_e, design)
    d_e = estimateGLMTagwiseDisp(d_e, design)
    fit1 = glmFit(d_e, design)
    res1 = residuals(fit1, type="deviance")
    ruv_cov = RUVr(round(d_e$counts), 
               as.character(rownames(d_e$counts)), 
               k=nruv, res1)
    uqobs = cbind(uqobs, ruv_cov$W)
    uqobs$Age <- NULL
    pathdf <- merge(pathdf, uqobs, by = "Sample", all.x = TRUE)
    
    rownames(pathdf) <-pathdf$barcode
    pathdf$barcode <- NULL
    pathdf = pathdf[colnames(mat),]
    print('RUV over')
    return (pathdf)
}
makeRegFormula = function(pathdf, path, nruv = NULL, int.var=NULL){
    flist = c('~', path)
    if (length(unique(pathdf$subtype_new)) > 1){
        flist = c(flist, '+ subtype_new')
    }
    if (nruv > 0){
        ruvw = paste0("W_", 1:nruv)
        flist = c(flist, " + ", paste(ruvw, collapse=" + "))
    }
    nb.form = asform(flist)
    print(nb.form)
    mdx = model.matrix(nb.form,data = pathdf)
    mdx = data.frame(mdx)
    leff = paste0('logFC_', 'Age')
    peff = paste0('p_', 'Age')
    return(list(nb.form=nb.form, mdx=mdx,
            leff=leff, peff=peff, pathstr='Age'))
}

makeNormMat = function(mat, csm){
    norm = as.matrix(mat)
    norm = sweep(norm, 2, csm / 10000,'/')
    norm = log(norm + 1)
    return(norm)
}


runSingleChunkMAST = function(mat, covmat, fdata, 
                              ind, csm = csm, nb.form, pathstr){
    # Make matrix and chunked sca object:
    norm = makeNormMat(mat[ind,], csm)
    sca = FromMatrix(exprsArray = norm, cData=covmat, fData=fdata[ind,,drop=F], check_sanity=TRUE)
    # Calculate regression, run LRT, refit:
    zlm.obj = zlm(nb.form, sca)
    # NOTE: if multiple pathstr, this will be multiple x slower
    summaryCond <- summary(zlm.obj, doLRT=pathstr) 
    # MAST results:
    summaryDt <- summaryCond$datatable
    resdf = NULL
    for (pstr in pathstr){
        fcHurdle <- merge(
            # Hurdle P values
            summaryDt[contrast==pstr & component=='H',
                .(primerid, `Pr(>Chisq)`)],
            # logFC coefficients
            summaryDt[contrast==pstr & component=='logFC',
                .(primerid, coef, ci.hi, ci.lo)], by='primerid')
        chunkdf = data.frame(fcHurdle)
        if (length(pathstr) == 1){
            colnames(chunkdf) = c('gene','p','coef','ci.hi','ci.lo')
        } else {
            colnames(chunkdf) = c('gene',
                paste0(c('p','coef','ci.hi','ci.lo'), '_', pstr))
        }
        if (is.null(resdf)){
            resdf = chunkdf
        } else {
            resdf = merge(resdf, chunkdf, all=TRUE)
        }
    }
    return(resdf)
}
region <- unique(ALL_rds$Region)
    regdir = "/hwfssz3/PS_JLU/zhangxiao6/Figure7/MBA_PFC_DEG_list"
  dir.create(paste(regdir, region, sep = '/'))
  dir.create(paste(regdir, region, 'RUV_result',sep = '/'))


 for (celltype in celltypes) {
        print(paste("Processing cell type:", celltype))
    
    # 1. 处理非法字符
     sanitized_celltype <- gsub("[^[:alnum:]]", "_", celltype)
    
    # 2. 生成文件名并检查是否已存在
    prefstr = paste(region, sanitized_celltype, sep = '_')
    file_name = paste0(prefstr,'_',method,'.csv')
    path_name = paste(regdir, region, file_name, sep = '/')
    
    # 如果文件已存在，跳过这个celltype
    if(file.exists(path_name)) {
        print(paste("文件已存在，跳过celltype:", celltype))
        next
    }
    
    rds <- subset(ALL_rds, ALL_rds$cluster1 == celltype)
    sample_counts <- table(rds@meta.data$Sample) 
    valid_samples <- names(sample_counts)[sample_counts >= 10]
    rds <- subset(rds, subset = Sample %in% valid_samples) 
    meta <- rds@meta.data
    meta_age_1 <- meta[meta$Age_group == 'Young',]
    meta_age_2 <- meta[meta$Age_group == 'Middle',]
    meta_age_3 <- meta[meta$Age_group == 'Old',]
    meta_age_4 <- meta[meta$Age_group == 'Very_old',]
    meta_list <- list(meta_age_1, meta_age_2, meta_age_3, meta_age_4)
    
    # 使用处理后的celltype名称
    prefstr = paste(region, sanitized_celltype, sep = '_')
    file_name = paste0(prefstr,'_',method,'.csv')
    ruv_name = paste0(prefstr,'_RUV.csv')
    ruv_path = paste(regdir, region, 'RUV_result', ruv_name, sep = '/')
    path_name = paste(regdir, region, file_name, sep = '/')
    outtsv = paste0(regdir, prefstr,'.tsv.gz')
    outrda = paste0(regdir, prefstr,'.rda')
    nsigfile = paste0(regdir, prefstr,'.nsig.tsv')
    
    print(method)
    pcut = 0.1
    print(prefstr)
    pathdf = rds@meta.data
    pathdf$barcode <- rownames(pathdf)
    mat = GetAssayData(rds, assay = "RNA", layer = "counts")
    path = 'Age'
    pctcells = rowSums(mat > 0) / ncol(mat)
    mat = subsetMatrixForDE(mat, pathdf, pctcells = pctcells, pctcut = pcut)
    nruv = 10
    uruv = 5
    pathdf = runRUVpsbulk(pathdf, mat, path, nruv = nruv)
    ll = makeRegFormula(pathdf, path=path, nruv = uruv)
    write.csv(pathdf, ruv_path)
    fdata = data.frame(primerid=rownames(mat))
    covmat = ll$mdx
    covmat$wellKey = rownames(pathdf)
    csm = colSums(mat)
    print("[STATUS] Running MAST")
    chunksize = 500
    nchunk = ceiling(nrow(mat) / chunksize)
    fulldf = c()
    t0 = proc.time()
    t1 = t0
    for (chunk in 1:nchunk){
        ind = (1 + (chunk-1) * chunksize):min(c(chunk * chunksize, nrow(mat))) 
        print(paste(chunk, length(ind)))
        chunkrda = paste0(regdir, prefstr, '.', 
                          chunksize, '_', chunk, '.rda')
        chunkdf = runSingleChunkMAST(mat=mat, covmat=covmat, fdata=fdata, 
                                         ind = ind, csm=csm, nb.form = ll$nb.form, pathstr=ll$pathstr)
        #save(chunkdf, file=chunkrda)  
        print(chunkrda)
        fulldf = rbind(fulldf, chunkdf)
        t1 = est.endtime(t0, t1, chunk, nchunk)
    }
    
    fulldf$fdr = p.adjust(fulldf$p, 'fdr')
    fulldf = fulldf[order(fulldf$p),]
    fulldf = fulldf[!is.na(fulldf$coef),]
    fulldf$subtype <- celltype
    fulldf$Region <- region
    fulldf$Stage <- 'pDEG'
    write.csv(fulldf, file = path_name)
}





