#' Identify over-expressed signaling genes associated with each cell group
#'
#' USERS can use customized gene set as over-expressed signaling genes by setting `object@var.features[[features.name]] <- features.sig`
#' The Bonferroni corrected/adjusted p value can be obtained via `object@var.features[[paste0(features.name, ".info")]]`. Note that by default `features.name = "features"`
#'
#' @param object CellChat object
#' @param data.use a customed data matrix. Default: data.use = NULL and the expression matrix in the slot 'data.signaling' is used
#' @param group.by cell group information; default is `object@idents`; otherwise it should be one of the column names of the meta slot
#' @param idents.use a subset of cell groups used for analysis
#' @param invert whether invert the idents.use
#' @param group.dataset dataset origin information in a merged CellChat object; set it as one of the column names of meta slot when identifying the highly enriched genes in one dataset for each cell group
#' @param pos.dataset the dataset name used for identifying highly enriched genes in this dataset for each cell group
#' @param features.name a char name used for storing the over-expressed signaling genes in `object@var.features[[features.name]]`
#' @param only.pos Only return positive markers
#' @param features features used for identifying Over Expressed genes. default use all features
#' @param return.object whether return the object; otherwise return a data frame consisting of over-expressed signaling genes associated with each cell group
#' @param thresh.pc Threshold of the percent of cells expressed in one cluster
#' @param thresh.fc Threshold of Log Fold Change
#' @param thresh.p Threshold of p-values
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom stats sd wilcox.test
#' @importFrom stats p.adjust
#'
#' @return A CellChat object or a data frame. If returning a CellChat object, two new elements named 'features.name' and paste0(features.name, ".info") will be added into the list `object@var.features`
#' `object@var.features[[features.name]]` is a vector consisting of the identified over-expressed signaling genes;
#' `object@var.features[[paste0(features.name, ".info")]]` is a data frame returned from the differential expression analysis
#' @export
#'
identifyOverExpressedGenes_chenge <- function(object, data.use = NULL, group.by = NULL, idents.use = NULL, invert = FALSE, group.dataset = NULL, pos.dataset = NULL, features.name = "features",  only.pos = TRUE, features = NULL, return.object = TRUE,
                                       thresh.pc = 0, thresh.fc = 0, thresh.p = 0.05) {
  if (!is.list(object@var.features)) {
    stop("Please update your CellChat object via `updateCellChat()`")
  }
  if (is.null(data.use)) {
    X <- object@data.signaling
    if (nrow(X) < 3) {stop("Please check `object@data.signaling` and ensure that you have run `subsetData` and that the data matrix `object@data.signaling` looks OK.")}
  } else {
    X <- data.use
  }

  if (is.null(features)) {
    features.use <- row.names(X)
  } else {
    features.use <- intersect(features, row.names(X))
  }
  data.use <- X[features.use,]
  data.use <- as.matrix(data.use)

  if (is.null(group.by)) {
    labels <- object@idents
    if (!is.factor(labels)) {
      message("Use the joint cell labels from the merged CellChat object")
      labels <- object@idents$joint
    }
  } else {
    labels <- object@meta[[group.by]]
  }
  if (!is.factor(labels)) {
    labels <- factor(labels)
  }
  level.use <- levels(labels)[levels(labels) %in% unique(labels)]
  if (!is.null(idents.use)) {
    if (invert) {
      level.use <- level.use[!(level.use %in% idents.use)]
    } else {
      level.use <- level.use[level.use %in% idents.use]
    }
  }
  numCluster <- length(level.use)

  if (!is.null(group.dataset)) {
    labels.dataset <- as.character(object@meta[[group.dataset]])
    if (!(pos.dataset %in% unique(labels.dataset))) {
      cat("Please set pos.dataset to be one of the following dataset names: ", unique(as.character(labels.dataset)))
      stop()
    }
  }

  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )

  mean.fxn <- function(x) {
    return(log(x = mean(x = expm1(x = x)) + 1))
  }
  labels <- as.character(labels)
  genes.de <- vector("list", length = numCluster)
  for (i in 1:numCluster) {
    features <- features.use
    if (is.null(group.dataset)) {
      cell.use1 <- which(labels == level.use[i]) #目标细胞类型
      cell.use2 <- base::setdiff(1:length(labels), cell.use1) #非目标的所有细胞
    } else {
      cell.use1 <- which((labels == level.use[i]) & (labels.dataset == pos.dataset))
      cell.use2 <- which((labels == level.use[i]) & (labels.dataset != pos.dataset))
    }

    # feature selection (based on percentages)
#    thresh.min <- 0
#    pct.1 <- round(
#      x = rowSums(data.use[features, cell.use1, drop = FALSE] > thresh.min) /
#        length(x = cell.use1),
#      digits = 3
#    )
#    pct.2 <- round(
#      x = rowSums(data.use[features, cell.use2, drop = FALSE] > thresh.min) /
#        length(x = cell.use2),
#      digits = 3
#    )
    thresh.min <- 0
    pct.1 <-  rowSums(data.use[features, cell.use1, drop = FALSE] > thresh.min) / length(x = cell.use1)
    pct.2 <- rowSums(data.use[features, cell.use2, drop = FALSE] > thresh.min) / length(x = cell.use2)
    
    data.alpha <- cbind(pct.1, pct.2)
    colnames(x = data.alpha) <- c("pct.1", "pct.2")
    alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
    names(x = alpha.min) <- rownames(x = data.alpha)
    features <- names(x = which(x = alpha.min > thresh.pc))
    if (length(x = features) == 0) {
      #stop("No features pass thresh.pc threshold")
      next
    }

    # feature selection (based on average difference)
    data.1 <- apply(X = data.use[features, cell.use1, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
    data.2 <- apply(X = data.use[features, cell.use2, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
    FC <- (data.1 - data.2) 
    if (only.pos) {
      features.diff <- names(which(FC > thresh.fc))  # 这里会去掉目标细胞类型里面平均量小于剩余细胞类型的平均表达量的gene
    } else {
      features.diff <- names(which(abs(FC) > thresh.fc))
    }

    #features <- intersect(x = features, y = features.diff)
    features <- unique(x = features, y = features.diff)
    if (length(x = features) == 0) {
      #  stop("No features pass thresh.fc threshold")
      next
    }

    data1 <- data.use[features, cell.use1, drop = FALSE]
    data2 <- data.use[features, cell.use2, drop = FALSE]

    pvalues <- unlist(
      x = my.sapply(
        X = 1:nrow(x = data1),
        FUN = function(x) {
          # return(wilcox.test(data1[x, ], data2[x, ], alternative = "greater")$p.value)
          return(wilcox.test(data1[x, ], data2[x, ])$p.value)
        }
      )
    )

    pval.adj = stats::p.adjust(
      p = pvalues,
      method = "bonferroni",
      n = nrow(X)
    )
    genes.de[[i]] <- data.frame(clusters = level.use[i], features = as.character(rownames(data1)), pvalues = pvalues, logFC = FC[features], data.alpha[features,, drop = F],pvalues.adj = pval.adj, stringsAsFactors = FALSE)
  }

  markers.all <- data.frame()
  for (i in 1:numCluster) {
    gde <- genes.de[[i]]
    if (!is.null(gde)) {
      gde <- gde[order(gde$pvalues, -gde$logFC), ]
      gde <- subset(gde, subset = pvalues < thresh.p)
      if (nrow(gde) > 0) {
        markers.all <- rbind(markers.all, gde)
      }
    }
  }
  if (only.pos & nrow(markers.all) > 0) {
    #markers.all <- subset(markers.all, subset = logFC > 0)
    markers.all <- markers.all
  }
  if (!is.null(group.dataset)) {
    markers.all$datasets[markers.all$logFC > 0] <- pos.dataset
    markers.all$datasets[markers.all$logFC < 0] <- setdiff(unique(labels.dataset), pos.dataset)
    markers.all$datasets <- factor(markers.all$datasets, levels = levels(factor(object@meta[[group.dataset]])))
    markers.all <- markers.all[order(markers.all$datasets, markers.all$pvalues, -markers.all$logFC), ]
  }
  markers.all$features <- as.character(markers.all$features)

  features.sig <- markers.all$features
  object@var.features[[features.name]] <- features.sig
  features.name <- paste0(features.name, ".info")
  object@var.features[[features.name]] <- markers.all

  if (return.object) {
    return(object)
  } else {
    return(markers.all)
  }
}
