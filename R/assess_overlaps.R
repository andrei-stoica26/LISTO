#' Assess the overlap of two lists of marker data frames.
#'
#' This function assesses the overlap of two lists of marker data frames.
#'
#' @param markerList1 List of data frames.
#' @param markerList2 List of data frames.
#' @inheritParams markerSetsPhyper
#' @param pvalThr p-value threshold.
#' @param verbose Whether the output should be verbose.
#'
#' @return A data frame.
#'
#' @export
#'
markerListPhyper <- function(markerList1, markerList2, nGenes,
                             colStr = 'avg_log2FC',
                             isHighTop = TRUE,
                             extraCutoff = 0,
                             pvalThr = 0.05,
                             verbose = TRUE){
    df <- pairDF(names(markerList1), names(markerList2))
    df$pval <- apply(df, 1,
                     function(x){
                         markerNames1 <- x[1]
                         markerNames2 <- x[2]
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

                         }
                     )
    df <- byCorrectDF(df, pvalThr)
    return(df)
}
