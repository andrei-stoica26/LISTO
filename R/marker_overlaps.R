#' @importFrom parallel clusterExport makeCluster parSapply stopCluster
#'
NULL

#' Generate cutoffs for assessing marker overlaps
#'
#' This function generates cutoffs for assessing marker overlaps.
#'
#' @param markers1 Data frame of markers with at least one numeric column.
#' @param markers2 Data frame of markers with at least one numeric column.
#' @param colStr Numeric column.
#' @param isHighTop Whether higher values in the numeric column correspond to
#' top markers.
#' @param extraCutoff Cutoff to be placed at one end of the cutoff list.
#' If \code{isHighTop} is \code{TRUE}, it must be lower than all the cutoffs
#' present in the marker list, and if \code{isHighTop} is \code{FALSE}, it must
#' be higher that all cutoffs present in the marker list.
#' @param maxNCutoffs Maximum number of cutoffs. If the two input data frames
#' contain more cutoffs than this value, only \code{maxNCutoffs} linearly
#' spaced cutoffs will be selected from the original cutoff list.
#' @param verbose Whether the output should be verbose.
#'
#' @return A numeric vector.
#'
#' @keywords internal
#'
generateCutoffs <- function(markers1,
                            markers2,
                            colStr = 'avg2_logFC',
                            isHighTop = TRUE,
                            extraCutoff = 0,
                            maxNCutoffs = 500,
                            verbose = FALSE){

    values1 <- markers1[[colStr]]
    values2 <- markers2[[colStr]]
    cutoffs <- unique(c(values1, values2))
    if (isHighTop)
        cutoffs <- cutoffs[cutoffs < min(max(values1), max(values2))] else
            cutoffs <- cutoffs[cutoffs > max(min(values1), min(values2))]
    cutoffs <- c(cutoffs, extraCutoff)
    cutoffs <- sort(cutoffs, decreasing=isHighTop)
    nCutoffs <- length(cutoffs)
    if (nCutoffs > maxNCutoffs){
        if(verbose)
            message('Too many cutoffs found in the input data frames. Only ',
                maxNCutoffs, ' will be used')
        cutoffs <- cutoffs[seq(1, nCutoffs, length.out=maxNCutoffs)]
    }
    return(cutoffs)
}

#' Compute the p-value of overlap for two marker data frames
#'
#' This function computes the p-value of overlap for two marker data frames.
#'
#' @inheritParams generateCutoffs
#' @param cutoffs Cutoffs for assessing marker overlaps.
#' @param nDatasets Number of datasets.
#' @param nCores Number of cores.
#' @param ... Additional parameters passed to \code{pvalOverlap} or
#' \code{pvalOverlapMN}.
#'
#' @return A p-value.
#'
#' @keywords internal
#'
markerDFPairPval <- function(markers1,
                             markers2,
                             colStr,
                             cutoffs,
                             nDatasets = c('1', '2'),
                             nCores = 1,
                             ...){

    if(nDatasets == 1)
        pvalFun <- eval(as.name('pvalOverlap')) else
            pvalFun <- eval(as.name('pvalOverlapMN'))

    if (nCores == 1){
        pvals <- vapply(cutoffs, function(cutoff){
            markerNames1 <- rownames(markers1[markers1[[colStr]] > cutoff, ])
            markerNames2 <- rownames(markers2[markers2[[colStr]] > cutoff, ])
            pval <- pvalFun(markerNames1, markerNames2, ...)
        }, numeric(1))

    } else{
        message('Computing in parallel. Using ', nCores, ' cores.')
        clust <- makeCluster(nCores)
        pvalFunArgs <- list(...)

        if (length(pvalFunArgs) == 1){
            allGenes <- pvalFunArgs[[1]]
            clusterExport(clust, c('markers1', 'markers2', 'cutoffs', 'allGenes'),
                          envir=environment())
            pvals <- parSapply(clust, cutoffs, function(cutoff){
                markerNames1 <- rownames(markers1[markers1[[colStr]] > cutoff, ])
                markerNames2 <- rownames(markers2[markers2[[colStr]] > cutoff, ])
                pval <- pvalFun(markerNames1, markerNames2, allGenes)
            })
        } else {
            allGenes1 <- pvalFunArgs[[1]]
            allGenes2 <- pvalFunArgs[[2]]
            clusterExport(clust, c('markers1', 'markers2', 'cutoffs',
                                   'allGenes1', 'allGenes2'), envir=environment())
            pvals <- parSapply(clust, cutoffs, function(cutoff){
                markerNames1 <- rownames(markers1[markers1[[colStr]] > cutoff, ])
                markerNames2 <- rownames(markers2[markers2[[colStr]] > cutoff, ])
                pval <- pvalFun(markerNames1, markerNames2, allGenes1, allGenes2)
            })
        }
        stopCluster(clust)
    }

    pval <- median(BY(pvals)$Adjusted.pvalues)
    return(pval)
}

#' Assess the overlap of two marker data frames.
#'
#' This function assesses the overlap of two marker data frames.
#'
#' @inheritParams generateCutoffs
#' @param genes1 Genes in the first dataset.
#' @param genes2 Genes in the second dataset. If \code{NULL} (as default), no
#' second dataset will be used.
#'
#' @return A numeric value (hypergeometric p-value).
#'
#' @export
#'
markerDFPairOverlap <- function(markers1,
                                markers2,
                                genes1,
                                genes2 = NULL,
                                nCores = 1,
                                colStr = 'avg_log2FC',
                                isHighTop = TRUE,
                                extraCutoff = 0,
                                maxNCutoffs = 500,
                                verbose = FALSE){

    cutoffs <- generateCutoffs(markers1, markers2, colStr, isHighTop,
                               extraCutoff, maxNCutoffs, verbose)

    if(is.null(genes2))
        return(markerDFPairPval(markers1,
                                markers2,
                                colStr,
                                cutoffs,
                                1,
                                nCores,
                                length(genes1)))

    return(markerDFPairPval(markers1,
                            markers2,
                            colStr,
                            cutoffs,
                            2,
                            nCores,
                            genes1,
                            genes2))

}

#' Filter marker list based on log2 fold-change and pct.1
#'
#' This function filters a marker list based on log2 fold-change and pct.1.
#'
#' @param markerList List of marker data frames.
#' @param logFCThr Value of average log2 fold-change above which markers will
#' be retained.
#' @param pct1Thr Value of fraction of cells expressing the marker above which
#' markers will be retained.
#'
#' @return A list of filtered marker data frames.
#'
#' @keywords internal
#'
filterMarkerList <- function(markerList, logFCThr = 0, pct1Thr = 0)
    return(lapply(filterMarkerList, function(x)
        x <- x[x$avg_log2FC > logFCThr & x$pct.1 > pct1Thr, ]))

#' Assess the overlap of two lists of marker data frames.
#'
#' This function assesses the overlap of two lists of marker data frames.
#'
#' @param markerList1 List of marker data frames.
#' @param markerList2 List of marker data frames.
#' @inheritParams markerDFPairOverlap
#' @inheritParams filterMarkerList
#' @param pvalThr p-value threshold to be used by the Bonferroni correction.
#'
#' @return A data frame.
#'
#' @export
#'
markerDFListOverlap <- function(markerList1,
                                markerList2,
                                genes1,
                                genes2 = NULL,
                                nCores = 1,
                                logFCThr = 0,
                                pct1Thr = 0,
                                colStr = 'avg_log2FC',
                                isHighTop = TRUE,
                                extraCutoff = 0,
                                maxNCutoffs = 500,
                                pvalThr = 0.05,
                                verbose = FALSE){

    df <- expand.grid(names(markerList1), names(markerList2))
    colnames(df) <- c('Group1', 'Group2')
    if (logFCThr > 0 | pct1Thr > 0){
        markerList1 <- filterMarkerList(markerList1, logFCThr, pct1Thr)
        markerList2 <- filterMarkerList(markerList2, logFCThr, pct1Thr)
    }

    df$pval <- vapply(seq(nrow(df)),
                      function(i){
                          markerNames1 <- df[i, 1]
                          markerNames2 <- df[i, 2]
                          message('Assessing overlap between marker sets: ',
                                markerNames1, ' and ', markerNames2, '...')
                          markerDFPairOverlap(
                              markerList1[[markerNames1]],
                              markerList2[[markerNames2]],
                              genes1,
                              genes2,
                              nCores,
                              colStr,
                              isHighTop,
                              extraCutoff,
                              maxNCutoffs,
                              verbose)
                          }, numeric(1))
    df <- byCorrectDF(df, pvalThr)
    return(df)
}
