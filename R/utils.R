#' @importFrom hammers bfCorrectDF byCorrectDF
#' @importFrom sgof BY
#' @importFrom stats median phyper
#'
NULL

#' Calculate the hypergeometric p-value of enrichment for two sets
#'
#' This function calculates the hypergeometric p-value of enrichment for
#' two sets
#'
#' @param a Character vector
#' @param b Chracter vector
#' @param n Number of elements of set from which \code{a} and \code{b} are
#' selected
#' @param lowerTail Whether to calculated underenrichment (\code{TRUE}) or
#' overenrichment (\code{FALSE})
#'
#' @export
#'
setsPhyper <- function(a, b, n, lowerTail=FALSE){
    na <- length(a)
    nb <- length(b)
    nShared <- length(intersect(a, b))
    return (phyper(nShared - (1 - lowerTail), na, n - na, nb,
                   lower.tail=lowerTail))
}

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
#'
#' @return A numeric vector.
#'
#' @keywords internal
#'
generateCutoffs <- function(markers1,
                            markers2,
                            colStr = 'avg2_logFC',
                            isHighTop = TRUE,
                            extraCutoff = 0){
    values1 <- markers1[[colStr]]
    values2 <- markers2[[colStr]]
    cutoffs <- unique(c(values1, values2))
    if (isHighTop)
        cutoffs <- cutoffs[cutoffs < min(max(values1), max(values2))] else
            cutoffs <- cutoffs[cutoffs > max(min(values1), min(values2))]
    cutoffs <- sort(c(cutoffs, extraCutoff), decreasing=isHighTop)
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
                             extraCutoff = 0){
    cutoffs <- generateCutoffs(markers1, markers2, colStr, isHighTop,
                               extraCutoff)
    pvals <- vapply(cutoffs, function(cutoff){
        markerNames1 <- rownames(markers1[markers1[[colStr]] > cutoff, ])
        markerNames2 <- rownames(markers2[markers2[[colStr]] > cutoff, ])
        pval <- setsPhyper(markerNames1, markerNames2, nGenes)
    }, numeric(1))
    pval <- median(BY(pvals)$Adjusted.pvalues)
    return(pval)
}

#' Create a data frame from vector pairs
#'
#' This function creates a data frame from vector pairs.
#'
#' @param v Character vector.
#' @param w Character vector.
#'
#' @return A data frame with two columns.
#'
#' @noRd
#'
pairDF <- function(v, w){
    df <- data.frame(Group1 = unlist(lapply(v, function(x) rep(x, length(w)))),
                     Group2 = rep(w, length(v)))
    return(df)
}

