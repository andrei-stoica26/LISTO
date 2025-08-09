#' Generate cutoffs for assessing marker overlaps
#'
#' This function generates cutoffs for assessing marker overlaps.
#'
#' @param markers1 Data frame.
#' @param markers2 Data frame.
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
                            maxNCutoffs = 10000){
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
        message('Too many cutoffs found in the input data frames. Only ',
                maxNCutoffs, ' will be used')
        cutoffs <- cutoffs[seq(1, nCutoffs, length.out=maxNCutoffs)]
    }
    return(cutoffs)
}

#' Assess the overlap of two marker data frames.
#'
#' This function assesses the overlap of two marker data frames.
#'
#' @inheritParams generateCutoffs
#' @param nGenes Number of genes in the dataset.
#'
#' @return A numeric value (hypergeometric p-value).
#'
#' @export
#'
markerSetsPhyper <- function(markers1, markers2, nGenes,
                             colStr = 'avg_log2FC',
                             isHighTop = TRUE,
                             extraCutoff = 0,
                             maxNCutoffs = 10000){
    cutoffs <- generateCutoffs(markers1, markers2, colStr, isHighTop,
                               extraCutoff, maxNCutoffs)
    pvals <- vapply(cutoffs, function(cutoff){
        markerNames1 <- rownames(markers1[markers1[[colStr]] > cutoff, ])
        markerNames2 <- rownames(markers2[markers2[[colStr]] > cutoff, ])
        pval <- setsPhyper(markerNames1, markerNames2, nGenes)
    }, numeric(1))
    pval <- median(BY(pvals)$Adjusted.pvalues)
    return(pval)
}


#' Assess the overlap of two lists of marker data frames.
#'
#' This function assesses the overlap of two lists of marker data frames.
#'
#' @param markerList1 List of data frames.
#' @param markerList2 List of data frames.
#' @inheritParams markerSetsPhyper
#' @param logFCThr Value of average log2 fold-change above which markers will
#' be retained.
#' @param pct1Thr Value of fraction of cells expressing the marker above which
#' markers will be retained.
#' @inheritParams filterMarkerList
#' @param pvalThr p-value threshold to be used by the Bonferroni correction.
#' @param verbose Whether the output should be verbose.
#'
#' @return A data frame.
#'
#' @export
#'
markerListPhyper <- function(markerList1, markerList2, nGenes,
                             logFCThr = 0,
                             pct1Thr = 0,
                             colStr = 'avg_log2FC',
                             isHighTop = TRUE,
                             extraCutoff = 0,
                             maxNCutoffs = 10000,
                             pvalThr = 0.05,
                             verbose = TRUE){

    df <- expand.grid(names(markerList1), names(markerList2))
    if (logFCThr > 0 | pct1Thr > 0){
        markerList1 <- filterMarkerList(markerList1, logFCThr, pct1Thr)
        markerList2 <- filterMarkerList(markerList2, logFCThr, pct1Thr)

        nLists <- length(markerList1)
        markerList1 <- lapply(nLists, function(i){
            x <- markerList1[[i]]
            x <- x[x$avg_log2FC > logFCThr & x$pct.1 > pct1Thr, ]
            if (verbose)
                message(nrow(x), ' markers retained for ',
                        names(markerList1)[1], '.')
            return(x)
        })
    }

    df$pval <- vapply(seq(nrow(df)),
                      function(i){
                          markerNames1 <- df[i, 1]
                          markerNames2 <- df[i, 2]
                          if (verbose)
                              message('Assessing overlap between marker sets: ',
                                     markerNames1, ' and ', markerNames2, '...')
                          markerSetsPhyper(
                              markerList1[[markerNames1]],
                              markerList2[[markerNames2]],
                              nGenes,
                              colStr,
                              isHighTop,
                              extraCutoff,
                              maxNCutoffs)
                          }, numeric(1))
    df <- byCorrectDF(df, pvalThr)
    return(df)
}
