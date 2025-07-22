#' Assess the overlap of two lists of marker data frames.
#'
#' This function assesses the overlap of two lists of marker data frames.
#'
#' @param markerList1 List of data frames.
#' @param markerList2 List of data frames.
#' @inheritParams markerSetsPhyper
#' @param pvalThr p-value threshold.
#'
#' @return A data frame.
#'
#' @export
#'
markerListPhyper <- function(markerList1, markerList2, nGenes,
                             colStr = 'avg_log2FC',
                             isHighTop = TRUE,
                             extraCutoff = 0,
                             pvalThr = 0.05){
    df <- pairDF(names(markerList1), names(markerList2))
    df$pval <- apply(df, 1,
                     function(x) markerSetsPhyper(
                         markerList1[[x[1]]],
                         markerList2[[x[2]]],
                         nGenes,
                         colStr,
                         isHighTop,
                         extraCutoff)
                     )
    df <- byCorrectDF(df, pvalThr)
    return(df)
}
