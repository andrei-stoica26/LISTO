#' Assess the overlap of two lists of marker data frames.
#'
#' This function assesses the overlap of two lists of marker data frames.
#'


#' @param logFCThr Value of average log2 fold-change above which markers will
#' be retained.
#' @param pct1Thr Value of fraction of cells expressing the marker above which
#' markers will be retained.



#' Assess the overlap of two lists of marker data frames.
#'
#' This function assesses the overlap of two lists of marker data frames.
#'
#' @param markerList1 List of data frames.
#' @param markerList2 List of data frames.
#' @inheritParams markerSetsPhyper
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
                              extraCutoff)
                          }, numeric(1))
    df <- byCorrectDF(df, pvalThr)
    return(df)
}
