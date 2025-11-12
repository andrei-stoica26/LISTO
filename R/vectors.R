#' Add numeric vectors of different lenghts
#'
#' This function adds numeric vectors of different lengths by filling shorter
#' vectors with zeroes.
#'
#' @param ... Numeric vectors.
#'
#' @return A numeric vector.
#'
#' @examples
#' vSum(c(1, 4), c(2, 8, 6), c(1, 7), c(10, 4, 6, 7))
#'
#' @export
#'
vSum <- function(...){
    vectors <- list(...)
    lengths <- vapply(vectors, length, integer(1))
    maxLen <- max(lengths)
    vectors <- mapply(function(v, l) c(v, rep(0, maxLen - l)), vectors,
                      lengths, SIMPLIFY=FALSE)
    return(psum(do.call(psum, vectors)))
}

#' Compute the prime factor decomposition of the binomial coefficient
#'
#' This function computes the prime factor decomposition of the
#' binomial coefficient.
#'
#' @param n Total number of elements.
#' @param k Number of selected elements.
#'
#' @return A vector in which positions represent prime numbers and values
#' represent their exponents in the factorial decomposition.
#'
#' @examples
#' choose(8, 4)
#'
#' @export
#'
vChoose <- function(n, k)
    return(vSum(factorialPrimePowers(n),
                -1 * factorialPrimePowers(k),
                -1 * factorialPrimePowers(n - k)))

#' Compute the prime representation of the numerator of the fraction
#' representing the probability that two subsets of sets M and N intersect
#' in k points
#'
#' This function computes the numerator of the fraction representing the
#' probability that two subsets of sets M and N intersect in k points
#'
#' @param intMN Number of elements in the intersection of sets M and N.
#' @param intAN Number of elements in the intersection of sets A (subset of M)
#' and N.
#' @param intBM Number of elements in the intersection of sets B (subset of N)
#' and M.
#' @param k Number of elements in the intersection of sets A and B.
#'
#' @return A vector containing the prime representation of the fraction
#' representing the probability that two subsets of sets M and N intersect in k
#' points. Positions represent prime numbers in order (2, 3, 5...), and values
#' represent their exponents in the prime decomposition.
#'
#' @keywords internal
#'
vNumeratorMN <- function(intMN, intAN, intBM, k)
    return(vSum(vChoose(intAN, k), vChoose(intMN - intAN, intBM - k)))
