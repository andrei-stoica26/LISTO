#' Compute the probability that two subsets of sets M and N intersect in
#' at least k points
#'
#' This function computes the probability that two subsets of sets M and N
#' intersect in at least k points.
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
pvalCountsMN <- function(intMN, intAN, intBM, k){
    denom <- -1 * vChoose(intMN, intBM)
    pval <- sum(vapply(seq(k, min(intAN, intBM)), function(i){
        exponents <- vSum(vNumeratorMN(intMN, intAN, intBM, i), denom)
        primes <- generate_n_primes(length(exponents))
        return(powerProduct(primes, exponents))
    }, numeric(1)))
    return(pval)
}
