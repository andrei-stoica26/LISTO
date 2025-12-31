#' @importFrom methods is
#' @importFrom parallel clusterExport makeCluster parSapply stopCluster
#'
NULL

#' Filter items based on a provided cutoff
#'
#' This function filters items based on a provided cutoff.
#'
#' @inheritParams getObjectValues
#' @param cutoff Cutoff for assessing item overlaps.
#' @param compFun Comparison function.
#'
#' @keywords internal
#'
filterItems <- function(obj, numCol = NULL, cutoff = NULL, compFun = `>`){
    if (is.null(numCol) | is(obj, 'character'))
        return(obj)
    return(rownames(obj[compFun(obj[[numCol]], cutoff), ]))
}

#' Compute the p-value of overlap for two marker data frames
#'
#' This function computes the p-value of overlap for two marker data frames.
#'
#' @inheritParams generateCutoffs
#' @param universe1 The set from which the items stored
#' in \code{obj1} are selected.
#' @param universe2 The set from which the items stored
#' in \code{obj2} are selected.
#'
#' @return A p-value.
#'
#' @keywords internal
#'
pvalObjectsCore <- function(obj1,
                            obj2,
                            obj3 = NULL,
                            universe1,
                            universe2 = NULL,
                            numCol = NULL,
                            cutoff = NULL,
                            compFun = `>`,
                            type = c('2N', '2MN', '3N')){

    names1 <- filterItems(obj1, numCol, cutoff, compFun)
    names2 <- filterItems(obj2, numCol, cutoff, compFun)
    if (type == '2N')
        return(pvalSets2N(names1, names2, universe1))

    if (type == '2MN')
        return(pvalSets2MN(names1, names2, universe1, universe2))

    if (type == '3N'){
        names3 <- filterItems(obj3, numCol, cutoff, compFun)
        return(pvalSets3N(names1, names2, names3, universe1))
    }
}

#' Assess the overlap of two marker data frames.
#'
#' This function assesses the overlap of two marker data frames.
#'
#' @inheritParams generateCutoffs
#' @inheritParams pvalObjectsCore
#' @param mtMethod Multiple testing correction method.
#' @param nCores Number of cores.
#' @param type Type of overlap assessment. Choose between: two sets belonging
#' to the same universe ('2N'), two sets belonging to different universes
#' ('2MN'), three sets belonging to the same universe ('3MN').
#'
#' @return A numeric value (p-value).
#'
#' @export
#'
pvalObjects <- function(obj1,
                        obj2,
                        obj3 = NULL,
                        universe1,
                        universe2 = NULL,
                        numCol = NULL,
                        isHighTop = TRUE,
                        maxCutoffs = 500,
                        mtMethod = c('BY', 'holm', 'hochberg',
                                     'hommel', 'bonferroni', 'BH',
                                     'fdr', 'none'),
                        nCores = 1,
                        type = c('2N', '2MN', '3N')){

    mtMethod <- match.arg(mtMethod, c('BY', 'holm', 'hochberg',
                                      'hommel', 'bonferroni', 'BH',
                                      'fdr', 'none'))
    type <- match.arg(type, c('2N', '2MN', '3N'))

    cutoffs <- generateCutoffs(obj1, obj2, obj3, numCol, isHighTop,
                               maxCutoffs)
    compFun <- ifelse(isHighTop, `>`, `<`)

    if (nCores == 1){
        pvals <- vapply(cutoffs, function(cutoff)
            return(pvalObjectsCore(obj1, obj2, obj3, universe1, universe2,
                                   numCol, cutoff, compFun, type)), numeric(1))

    } else{
        clust <- makeCluster(nCores)
        clusterExport(clust,
                      c('obj1', 'obj2', 'obj3',  'universe1', 'universe2',
                        'numCol', 'cutoffs', 'compFun', 'type'),
                      envir=environment())
        pvals <- parSapply(clust, cutoffs, function(cutoff)
            return(pvalObjectsCore(obj1, obj2, obj3, universe1, universe2,
                                   numCol, cutoff, compFun, type)))
        stopCluster(clust)
    }

    pval <- mtCorrectV(pvals, mtMethod, 'median')
    return(pval)
}
