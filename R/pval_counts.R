#' Compute the probability that two subsets of sets M and N intersect in
#' at least k points
#'
#' This function computes the probability that two subsets A and B of sets
#' M and N intersect in at least k points.
#'
#' @inheritParams vNumeratorMN
#'
#' @return A p-value.
#'
#' @examples
#' pvalCounts2MN (300, 23, 24, 6)
#'
#' @export
#'
pvalCounts2MN <- function(intMN, intAN, intBM, k){
    denom <- -1 * vChoose(intMN, intBM)
    pval <- sum(vapply(seq(k, min(intAN, intBM)), function(i){
        exponents <- vSum(vNumeratorMN(intMN, intAN, intBM, i), denom)
        primes <- generate_n_primes(length(exponents))
        return(powerProduct(primes, exponents))
    }, numeric(1)))
    return(pval)
}

#' Compute the probability that three subsets of a set intersect in at least k points
#'
#' This function computes the probability that three subsets of a set intersect
#' in at least k points.
#'
#' @param lenA Size of the first subset.
#' @param lenB Size of the second subset.
#' @param lenC Size of the third subset.
#' @param n Size of the set comprising the subsets.
#' @param k Size of the intersection.
#'
#' @return A p-value.
#'
#' @examples
#' pvalCounts3N (300, 200, 250, 400, 180)
#'
#' @export
#'
pvalCounts3N <- function(lenA, lenB, lenC, n, k){
    if (k == 0)
        return(1)
    lengths <- sort(c(lenA, lenB, lenC))
    if (k > lengths[[1]])
        stop('`k` must not exceed the minimum length of the three sets.')
    return (sum(vapply(seq(k, lengths[[1]]),
                       function(x) probCounts3N(lengths[[1]],
                                                lengths[[2]],
                                                lengths[[3]],
                                                n,
                                                x), numeric(1))))
}

