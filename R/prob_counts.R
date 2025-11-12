#' @importFrom primes generate_n_primes generate_primes
#' @importFrom statisfactory psum
#'
NULL

#' Compute the probability that two subsets of sets M and N intersect in k
#' points
#'
#' This function computes the probability that two subsets of sets M and N
#' intersect in k points.
#'
#' @inheritParams vNumeratorMN
#'
#' @return The probability that two subsets of sets M and N intersect in k
#' points.
#'
#' @examples
#' probCountsMN(8, 6, 4, 2)
#'
#' @export
#'
probCountsMN <- function(intMN, intAN, intBM, k){
    exponents <- vSum(vNumeratorMN(intMN, intAN, intBM, k),
                      -1 * vChoose(intMN, intBM))
    primes <- generate_n_primes(length(exponents))
    return(powerProduct(primes, exponents))
}
