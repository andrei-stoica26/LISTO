#' @importFrom primes generate_n_primes generate_primes
#' @importFrom statisfactory psum
#'
NULL

#' Compute the probability that two subsets of sets M and N intersect in k
#' points
#'
#' This function computes the probability that two subsets of sets M and N
#' intersect in k points. Intersection sizes (M with N, A with N and B with M)
#' must be provided.
#'
#' @inheritParams vNumeratorMN
#'
#' @return A numeric value in [0, 1].
#'
#' @examples
#' probCounts2MN(8, 6, 4, 2)
#'
#' @export
#'
probCounts2MN <- function(intMN, intAN, intBM, k){
    exponents <- vSum(vNumeratorMN(intMN, intAN, intBM, k),
                      -1 * vChoose(intMN, intBM))
    primes <- generate_n_primes(length(exponents))
    return(powerProduct(primes, exponents))
}

#' Compute the probability that three subsets of given sizes intersect in k
#' points
#'
#' This function computes the probability that three subsets of of given sizes
#' intersect in k points.
#'
#' @param a Size of the first subset.
#' @param b Size of the second subset.
#' @param c Size of the third subset.
#' @param n Size of the set.
#' @param k Size of the intersection.
#'
#' @return A numeric value in [0, 1].
#'
#' @examples
#' probCounts3N(8, 6, 10, 20, 3)
#'
#' @export
#'
probCounts3N <- function(a, b, c, n, k){
    v <- seq(max(a + b - n, k), min(a, b, n + k - c))
    prob <- sum(vapply(v, function(x) dhyper(x, a, n - a, b) *
                            dhyper(k, x, n - x, c), numeric(1)))
    return (prob)
}
