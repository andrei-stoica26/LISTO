#' @importFrom stats median p.adjust phyper
#'
NULL

#' Calculate the hypergeometric p-value of enrichment for two sets
#'
#' This function calculates the hypergeometric p-value of enrichment for
#' two sets.
#'
#' @details This function calculates the hypergeometric p-value of enrichment
#' for two sets.
#'
#' @param a Character vector.
#' @param b Chracter vector.
#' @param n Number of elements of set from which \code{a} and \code{b} are
#' selected.
#' @param lowerTail Whether to calculate underenrichment (\code{TRUE}) or
#' overenrichment (\code{FALSE}).
#'
#' @export
#'
pvalOverlap <- function(a, b, n, lowerTail=FALSE){
    na <- length(a)
    nb <- length(b)
    nShared <- length(intersect(a, b))
    return (phyper(nShared - (1 - lowerTail),
                   na,
                   n - na,
                   nb,
                   lower.tail=lowerTail))
}

#' Compute the p-value of intersection of two gene subsets of sets M and N
#'
#' This function computes the p-value of intersection of two gene subsets
#' of sets M and N.
#'
#' @details A thin wrapper around \code{pvalOverlapMNk}.
#'
#' @inheritParams pOverlapMN
#'
#' @return The probability of intersection of the two gene subsets.
#'
#' @examples
#' pvalOverlapMN(LETTERS[seq(4, 10)],
#' LETTERS[seq(7, 15)],
#' LETTERS[seq(19)],
#' LETTERS[seq(6, 26)])
#'
#' @export
#'
pvalOverlapMN <- function(a, b, m, n){
    if (length(setdiff(a, m)))
        stop('`a` must be a subset of `m`.')
    if (length(setdiff(b, n)))
        stop('`b` must be a subset of `n`.')
    return(pvalOverlapMNk(length(intersect(m, n)),
                          length(intersect(a, n)),
                          length(intersect(b, m)),
                          length(intersect(a, b))
    ))
}
