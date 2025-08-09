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



