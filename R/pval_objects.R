#' @importFrom parallel clusterExport makeCluster parSapply stopCluster
#'
NULL


#' Compute the p-value of overlap for two marker data frames
#'
#' This function computes the p-value of overlap for two marker data frames.
#'
#' @inheritParams generateCutoffs
#' @param cutoffs Cutoffs for assessing item overlaps.
#' @param nDatasets Number of datasets. Choose between 1 (the two objects are
#' selected from the same dataset) and 2 (the two objects are selected
#' from different datasets).
#' @param nCores Number of cores. Only used if \code{nDatasets} is 2.
#' @param mtMethod Multiple testing correction method.
#' @param ... Additional parameters passed to \code{pvalOverlap} or
#' \code{pvalOverlapMN}.
#'
#' @return A p-value.
#'
#' @keywords internal
#'
pvalObjectsCore <- function(obj1,
                            obj2,
                            col,
                            cutoffs,
                            nDatasets = 1,
                            nCores = 1,
                            mtMethod = c('BY', 'holm', 'hochberg',
                                         'hommel', 'bonferroni', 'BH',
                                         'fdr', 'none'),
                            ...){

    mtMethod <- match.arg(mtMethod, c('BY', 'holm', 'hochberg',
                                      'hommel', 'bonferroni', 'BH',
                                      'fdr', 'none'))

    if(nDatasets %in% c(1, 2))
        stop('`nDatasets` must be either 1 or 2.')

    if(nDatasets == 1)
        pvalFun <- eval(as.name('pvalSubsetsN')) else
            pvalFun <- eval(as.name('pvalSubsetsMN'))

        if (nCores == 1 | nDatasets == 1){
            pvals <- vapply(cutoffs, function(cutoff){
                names1 <- rownames(obj1[obj1[[col]] > cutoff, ])
                names2 <- rownames(obj2[obj2[[col]] > cutoff, ])
                pval <- pvalFun(names1, names2, ...)
            }, numeric(1))

        } else{
            clust <- makeCluster(nCores)
            pvalFunArgs <- list(...)
            allItems1 <- pvalFunArgs[[1]]
            allItems2 <- pvalFunArgs[[2]]
            clusterExport(clust, c('obj1', 'obj2', 'cutoffs',
                                   'allItems1', 'allItems2'), envir=environment())
            pvals <- parSapply(clust, cutoffs, function(cutoff){
                names1 <- rownames(obj1[obj1[[col]] > cutoff, ])
                names2 <- rownames(obj2[obj2[[col]] > cutoff, ])
                pval <- pvalFun(names1, names2, allItems1, allItems2)
            })
            stopCluster(clust)
        }

        pval <- mtCorrectV(pvals, mtMethod, 'median')
        return(pval)
}

#' Assess the overlap of two marker data frames.
#'
#' This function assesses the overlap of two marker data frames.
#'
#' @inheritParams generateCutoffs
#' @inheritParams pvalObjectsCore
#' @param allItems1 Either items in the first dataset or their number. If a
#' second dataset is provided (that is, \code{allItems} is not \code{NULL}),
#' the items in the first dataset rather than their number must be provided
#' here.
#' @param allItems2 All items in the second dataset. If \code{NULL}
#' (as default), no second dataset will be used.
#'
#' @return A numeric value (hypergeometric p-value).
#'
#' @export
#'
pvalObjects <- function(obj1,
                        obj2,
                        col,
                        allItems1,
                        allItems2 = NULL,
                        nCores = 1,
                        isHighTop = TRUE,
                        extraCutoff = 0,
                        maxNCutoffs = 500,
                        mtMethod = c('BY', 'holm', 'hochberg',
                                     'hommel', 'bonferroni', 'BH',
                                     'fdr', 'none'),
                        verbose = FALSE){

    mtMethod <- match.arg(mtMethod, c('BY', 'holm', 'hochberg',
                                      'hommel', 'bonferroni', 'BH',
                                      'fdr', 'none'))

    cutoffs <- generateCutoffs(obj1, obj2, col, isHighTop,
                               extraCutoff, maxNCutoffs, verbose)

    if(is.null(allItems2)){
        if(!is.numeric(allItems1))
            allItems1 <- length(allItems1)
        return(pvalObjectsCore(obj1, obj2, col, cutoffs, 1,
                               1, mtMethod, allItems1))
    }

    return(pvalObjectsCore(obj1, obj2, col, cutoffs, 2, nCores,
                           mtMethod, allItems1, allItems2))
}
