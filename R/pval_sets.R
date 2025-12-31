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
#' @param n Set from which \code{a} and \code{b} are
#' selected.
#'
#' @export
#'
pvalSets2N <- function(a, b, n){
    na <- length(a)
    nb <- length(b)
    nShared <- length(intersect(a, b))
    nn <- length(n)
    return (phyper(nShared - 1),
                   na,
                   nn - na,
                   nb,
                   lower.tail=FALSE)
}

#' Compute the p-value of intersection of two subsets of sets M and N
#'
#' This function computes the p-value of intersection of two subsets
#' of sets M and N.
#'
#' @details A thin wrapper around \code{pvalCounts2MN}.
#'
#' @inheritParams pvalSets2N
#' @param m Set from which \code{a} is selected.
#' @param n Set from which \code{b} is selected.
#'
#' @return A numeric value in [0, 1].
#'
#' @examples
#' pvalSets2MN(LETTERS[seq(4, 10)],
#' LETTERS[seq(7, 15)],
#' LETTERS[seq(19)],
#' LETTERS[seq(6, 26)])
#'
#' @export
#'
pvalSets2MN <- function(a, b, m, n){
    if (length(setdiff(a, m)))
        stop('`a` must be a subset of `m`.')
    if (length(setdiff(b, n)))
        stop('`b` must be a subset of `n`.')
    return(pvalCounts2MN(length(intersect(m, n)),
                         length(intersect(a, n)),
                         length(intersect(b, m)),
                         length(intersect(a, b))
    ))
}

#' Compute the p-value of intersection of three subsets
#'
#' This function computes the p-value of intersection of three subsets.
#'
#' @details A thin wrapper around \code{pvalCounts3N}.
#'
#' @inheritParams pvalSets2N
#' @param c A character vector.
#' @param n Set from which \code{a}, \code{b} and \code{c} are
#' selected.
#'
#' @return A numeric value in [0, 1].
#'
#' @examples
#' pvalSets3N(LETTERS[seq(4, 10)],
#' LETTERS[seq(7, 15)],
#' LETTERS[seq(19)],
#' LETTERS[seq(1, 26)])
#'
#' @export
#'
pvalSets3N <- function(a, b, c, n){
    if (length(setdiff(a, n)))
        stop('`a` must be a subset of `n`.')
    if (length(setdiff(b, n)))
        stop('`b` must be a subset of `n`.')
    if (length(setdiff(c, n)))
        stop('`c` must be a subset of `n`.')
    return(pvalCounts3N(length(a),
                        length(b),
                        length(c),
                        length(n),
                        length(Reduce(intersect, list(a, b, c)))))
}

