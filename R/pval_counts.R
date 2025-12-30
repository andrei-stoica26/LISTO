#' Compute the probability that two subsets of sets M and N intersect in
#' at least k points
#'
#' This function computes the probability that two subsets A and B of sets
#' M and N intersect in at least k points. In addition to k, the intersection
#' sizes of M and M, A and N, and B and M must pe provided as input.
#'
#' @inheritParams vNumeratorMN
#'
#' @return The probability that two subsets of sets M and N intersect in
#' at least k points.
#'
#' @examples
#' pvalCountsMN (300, 23, 24, 6)
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



TestThreeSetsExact <- function(a, b, c, N)
    #Must return 1
    return(sum(sapply(0:min(a,b,c), function(x) ThreeSetsExact(a, b, c, x, N))))
