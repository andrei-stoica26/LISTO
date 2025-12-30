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
#' @param a A character vector.
#' @param b A character vector.
#' @param n Number of elements of set from which \code{a} and \code{b} are
#' selected.
#'
#' @export
#'
pvalSubsetsN <- function(a, b, n){
    na <- length(a)
    nb <- length(b)
    nShared <- length(intersect(a, b))
    return (phyper(nShared - 1),
                   na,
                   n - na,
                   nb,
                   lower.tail=FALSE)
}

#' Compute the p-value of intersection of two subsets of sets M and N
#'
#' This function computes the p-value of intersection of two subsets
#' of sets M and N.
#'
#' @details A thin wrapper around \code{pvalOverlapMNk}.
#'
#' @inheritParams pvalSubsetsN
#' @param m Set from which \code{a} is selected.
#' @param n Set from which \code{b} is selected.
#'
#' @return The probability of intersection of the two subsets.
#'
#' @examples
#' pvalSubsetsMN(LETTERS[seq(4, 10)],
#' LETTERS[seq(7, 15)],
#' LETTERS[seq(19)],
#' LETTERS[seq(6, 26)])
#'
#' @export
#'
pvalSubsetsMN <- function(a, b, m, n){
    if (length(setdiff(a, m)))
        stop('`a` must be a subset of `m`.')
    if (length(setdiff(b, n)))
        stop('`b` must be a subset of `n`.')
    return(pvalCountsMN(length(intersect(m, n)),
                        length(intersect(a, n)),
                        length(intersect(b, m)),
                        length(intersect(a, b))
    ))
}

pvalThreeSubsetsN <- function(a, b, c, n){
    na <- length(a)
    nb <- length(b)
    nShared <- length(intersect(a, b))
    return (phyper(nShared - 1),
                   na,
                   n - na,
                   nb,
                   lower.tail=FALSE)
}
