library(Seurat)
library(MAST)
library(dplyr)
library(tidyr)
library(Matrix)
library(RUVSeq)
library(edgeR)
library(parallel)
library(dplyr)

# Load data
rds <- readRDS('~/Mic_snRNA_merge_annotation.rds')
rds

# Identify which column contains cell names
colnames(rds@meta.data)
unique(rds$cluster1)
rds$subtype <- rds$cluster1

# Define label list
label_list <- c('Mic1_PFC', 'Mic1_Hf','Mic1_Str','Mic1_Tha','Mic1_Hyt','Mic1_MB','Mic1_PM')
length(label_list)
rds@meta.data$label <- paste(rds@meta.data$subtype, rds@meta.data$Region, sep = '_')
length(intersect(unique(rds@meta.data$label), label_list))
label_list <- intersect(unique(rds@meta.data$label), label_list)
unique(rds@meta.data$label)
unique(rds@meta.data$Age_group)

# Helper function: create indicator matrix
make.tform = function(x, norm = FALSE, u = NULL) {
    if (is.null(u)){u = unique(x)}
    tmat = sapply(u, function(y){1 * (x == y)})
    if (norm){
        tmat = sweep(tmat, 2, apply(tmat, 2, sum), '/')
    }
    return(tmat)
}

# Helper function: create formula
asform = function(x){ as.formula(paste0(x, collapse = '')) }

# Helper function: estimate remaining time
est.endtime <- function(start, last, done, total) {
    curr = proc.time()
    last.step = (proc.time() - last)[3]
    elapsed = (proc.time() - start)[3]
    cat(paste0(done, '/', total, '\t'))
    cat(paste0(round(last.step, 1), 's\t'))
    cat(paste0(round(elapsed, 1), 's\t'))
    each.time = elapsed / done
    est.left = (total - done) * each.time
    cat(paste0('Left: ', round(est.left, 1), 's\t'))
    fin.time = Sys.time() + est.left
    cat(paste0('Est. ', format(fin.time, "%H:%M:%S"), '\n'))
    return(curr)
}

# Function: subset matrix for differential expression analysis
subsetMatrixForDE = function(mat, pathdf, pctcells, pctcut = NULL, selgenes = NULL) {
    if (is.null(selgenes)) {
        if (!is.null(pctcut)) {
            keep.genes = names(pctcells)[pctcells > pctcut]
        }
    } else {
        keep.genes = selgenes[selgenes %in% rownames(mat)]
    }
    mat = mat[keep.genes, pathdf$barcode]
    print(paste("[STATUS] Subsetting matrix to", 
                paste0(dim(mat), collapse = ' x '), '(genes x cells)'))
    return(mat)
}

# Function: run RUV on pseudobulk data
runRUVpsbulk = function(pathdf, mat, path, nruv) { 
    pids = as.character(unique(pathdf$Sample))
    tform = make.tform(pathdf$Sample, u = pids)
    data_ind = mat %*% tform 
    indvar = 'Sample'
    uqcols = c(indvar, path)
    dform = asform(c('~', path))
    uqobs = unique(pathdf[, uqcols])
    rownames(uqobs) = uqobs[[indvar]]
    uqobs = uqobs[colnames(data_ind),]
    design = model.matrix(dform, data = uqobs)
    d_e = DGEList(data_ind, genes = rownames(data_ind))
    keep = rowSums(cpm(d_e) > 1) >= 3
    d_e = d_e[keep, , keep.lib.sizes = FALSE]
    d_e = calcNormFactors(d_e, method = "TMM")
    d_e = estimateGLMCommonDisp(d_e, design)
    d_e = estimateGLMTagwiseDisp(d_e, design)
    fit1 = glmFit(d_e, design)
    res1 = residuals(fit1, type = "deviance")
    ruv_cov = RUVr(round(d_e$counts), 
                   as.character(rownames(d_e$counts)), 
                   k = nruv, res1)
    uqobs = cbind(uqobs, ruv_cov$W)
    # Remove stim column if exists
    if ("stim" %in% colnames(uqobs)) uqobs$stim <- NULL
    pathdf <- merge(pathdf, uqobs, by = "Sample", all.x = TRUE)
    
    rownames(pathdf) <- pathdf$barcode
    pathdf$barcode <- NULL
    pathdf = pathdf[colnames(mat),]
    print('RUV completed')
    return(pathdf)
}

# Function: create regression formula
makeRegFormula = function(pathdf, path, nruv = NULL, int.var = NULL) {
    flist = c('~', path)
    if (nruv > 0) {
        ruvw = paste0("W_", 1:nruv)
        flist = c(flist, " + ", paste(ruvw, collapse = " + "))
    }
    nb.form = asform(flist)
    print(nb.form)
    mdx = model.matrix(nb.form, data = pathdf)
    mdx = data.frame(mdx)
    leff = paste0('logFC_', path)
    peff = paste0('p_', path)
    return(list(nb.form = nb.form, mdx = mdx,
                leff = leff, peff = peff, pathstr = path))
}

# Function: create normalized matrix
makeNormMat = function(mat, csm) {
    norm = as.matrix(mat)
    norm = sweep(norm, 2, csm / 10000, '/')
    norm = log(norm + 1)
    return(norm)
}

# Function: run MAST analysis on a single chunk
runSingleChunkMAST = function(mat, covmat, fdata, 
                              ind, csm = csm, nb.form, pathstr) {
    norm = makeNormMat(mat[ind,], csm)
    sca = FromMatrix(exprsArray = norm, cData = covmat, fData = fdata[ind, , drop = F], check_sanity = TRUE)
    zlm.obj = zlm(nb.form, sca)
    summaryCond <- summary(zlm.obj, doLRT = pathstr) 
    summaryDt <- summaryCond$datatable
    resdf = NULL
    for (pstr in pathstr) {
        fcHurdle <- merge(
            summaryDt[contrast == pstr & component == 'H',
                      .(primerid, `Pr(>Chisq)`)],
            summaryDt[contrast == pstr & component == 'logFC',
                      .(primerid, coef, ci.hi, ci.lo)], by = 'primerid')
        chunkdf = data.frame(fcHurdle)
        if (length(pathstr) == 1) {
            colnames(chunkdf) = c('gene', 'p', 'coef', 'ci.hi', 'ci.lo')
        } else {
            colnames(chunkdf) = c('gene',
                                   paste0(c('p', 'coef', 'ci.hi', 'ci.lo'), '_', pstr))
        }
        if (is.null(resdf)) {
            resdf = chunkdf
        } else {
            resdf = merge(resdf, chunkdf, all = TRUE)
        }
    }
    return(resdf)
}

# Main analysis steps
regdir <- '~/Revise_DEG_LIST/new_Mic1_deg_list'
ruv_dir <- '~/Revise_DEG_LIST/new_Mic1_deg_list/Mic1_ruvlist'

# Set log file path
log_file <- paste0(regdir, "/analysis_log_", format(Sys.Date(), "%Y%m%d"), ".txt")

# Start logging (append mode, also output to console)
sink(file = log_file, append = TRUE, split = TRUE)

DefaultAssay(rds) <- 'RNA'
celltypes <- label_list
print("Cell type list:")
print(celltypes)
comparisons <- list( 
    "Old_vs_Exceptionally_old" = c("Old", "Exceptionally old")
)
method = 'RUV_MAST'
#region = unique(rds$Region)
for (celltype in celltypes) {
    print(paste("Processing cell type:", celltype))
    
    # Extract specific cell type
    cell_idx <- which(rds$label == celltype)
    celltype_rds <- subset(rds, cells = cell_idx)
    
    # Filter samples: at least 3 cells per sample
    sample_counts <- table(celltype_rds@meta.data$Sample)
    valid_samples <- names(sample_counts)[sample_counts >= 10]
    celltype_rds <- subset(celltype_rds, Sample %in% valid_samples)
    
    # Check if there are enough samples
    if (length(unique(celltype_rds$Sample)) < 10) {
        print(paste("Skipping", celltype, "- insufficient samples"))
        next
    }
    
    # Analyze each comparison group
    for (comp_name in names(comparisons)) {
        comp_groups <- comparisons[[comp_name]]
        
        # Check if comparison groups exist
        if (!all(comp_groups %in% unique(celltype_rds$Age_group))) {
            print(paste("Skipping", comp_name, "- groups not found"))
            next
        }
        
        # Extract comparison group data
        comp_rds <- subset(celltype_rds, Age_group %in% comp_groups)
        
        # Check if each group has enough samples
        group_sample_n <- sapply(comp_groups, function(g) {
            length(unique(comp_rds$Sample[comp_rds$Age_group == g]))
        })
        names(group_sample_n) <- comp_groups
        # Skip if any group has <2 samples
        if (any(group_sample_n < 2)) {
            print(paste("Skipping", comp_name, "- groups have <2 biological replicates:",
                        paste(names(group_sample_n)[group_sample_n < 2], collapse = ", ")))
            next
        }
        
        # Set comparison variable
        comp_rds$stim = 0
        comp_rds$stim[comp_rds$Age_group == comp_groups[2]] = 1
        
        # Prepare output filenames
        safe_celltype <- gsub("[^[:alnum:]]", "_", celltype)
        region <- unique(celltype_rds$Region)
        prefstr = paste(region, safe_celltype, comp_name, sep = '_')
        file_name = paste0(prefstr, '_', method, '_sDEG.csv')
        ruv_name = paste0(prefstr, '_RUV.csv')
        path_name = paste(regdir, file_name, sep = '/')
        ruv_path = paste(ruv_dir, ruv_name, sep = '/')
        
        # Skip if results already exist
        if (file.exists(path_name)) {
            print(paste("Skipping", comp_name, "- results already exist"))
            next
        }
        
        print(paste("Processing comparison:", comp_name))
        
        # Prepare data
        pathdf = comp_rds@meta.data
        pathdf$barcode <- rownames(pathdf)
        mat = GetAssayData(comp_rds, assay = "RNA", layer = "counts")
        
        # Filter lowly expressed genes
        pctcells = rowSums(mat > 0) / ncol(mat)
        pcut = 0.1
        mat = subsetMatrixForDE(mat, pathdf, pctcells = pctcells, pctcut = pcut)
        
        # RUV correction
        nruv = 10
        uruv = 5
        pathdf = runRUVpsbulk(pathdf, mat, path = "stim", nruv = nruv)
        
        w_columns <- grep("^W_", colnames(pathdf), value = TRUE)
        actual_nruv <- length(w_columns)
        
        # Set actual number of RUV factors to use
        uruv <- min(5, actual_nruv)  # Use at most 5 factors, but not more than available
        
        # Save RUV results
        write.csv(pathdf, ruv_path)
        
        # Prepare MAST analysis
        ll = makeRegFormula(pathdf, path = "stim", nruv = uruv)
        fdata = data.frame(primerid = rownames(mat))
        covmat = ll$mdx
        covmat$wellKey = rownames(pathdf)
        csm = colSums(mat)
        
        # Run MAST analysis
        print("[STATUS] Running MAST")
        chunksize = 500
        nchunk = ceiling(nrow(mat) / chunksize)
        fulldf = c()
        t0 = proc.time()
        t1 = t0
        
        for (chunk in 1:nchunk) {
            ind = (1 + (chunk - 1) * chunksize):min(c(chunk * chunksize, nrow(mat))) 
            print(paste("Chunk", chunk, "of", nchunk, "size:", length(ind)))
            
            chunkdf = runSingleChunkMAST(
                mat = mat, 
                covmat = covmat, 
                fdata = fdata, 
                ind = ind, 
                csm = csm, 
                nb.form = ll$nb.form, 
                pathstr = ll$pathstr
            )
            
            fulldf = rbind(fulldf, chunkdf)
            t1 = est.endtime(t0, t1, chunk, nchunk)
        }
        
        # Process results
        fulldf$fdr = p.adjust(fulldf$p, 'fdr')
        fulldf = fulldf[order(fulldf$p),]
        fulldf = fulldf[!is.na(fulldf$coef),]
        
        # Save results
        write.csv(fulldf, file = path_name)
        print(paste(Sys.time(), "Completed", comp_name, "for", celltype))
    }
}
sink()
