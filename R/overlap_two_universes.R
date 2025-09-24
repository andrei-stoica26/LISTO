#' @importFrom primes generate_n_primes generate_primes
#' @importFrom statisfactory psum
#'
NULL

#' Generate the prime factor decomposition of n factorial.
#'
#' This function generates the prime factor decomposition of n!
#'
#' @param n A positive integer.
#'
#' @return A vector in which positions represent prime numbers and values
#' represent their exponents in the factorial decomposition.
#'
#' @examples
#' factorialPrimePowers(8)
#'
#' @export
#'
factorialPrimePowers <- function(n){
    if (n %in% c(0, 1))
        return(NULL)
    if(n < 0)
        stop('`n` must be a non-negative integer.')
    primes <- generate_primes(max=n)
    nPrimes <- length(primes)
    result <- rep(0, nPrimes)
    for (i in seq(nPrimes)){
        k <- primes[i]
        while (k <= n){
            result[i] <- result[i] + floor(n / k)
            k <- k * primes[i]
        }
    }
    return(result)
}

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
#' powerProduct(c(2, 3, 5), c(4, 2, 6))
#'
#' @export
#'
powerProduct <- function(primes, exponents)
    return(prod(mapply(function(x, y) x ^ y, primes, exponents)))


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
vSubsetOverlapProbMNkNumerator <- function(intMN, intAN, intBM, k)
    return(vSum(vChoose(intAN, k), vChoose(intMN - intAN, intBM - k)))

#' Compute the probability that two subsets of sets M and N intersect in k
#' points
#'
#' This function computes the probability that two subsets of sets M and N
#' intersect in k points.
#'
#' @inheritParams vSubsetOverlapProbMNkNumerator
#'
#' @return The probability that two subsets of sets M and N intersect in k
#' points.
#'
#' @examples
#' subsetOverlapProbMNk(8, 6, 4, 2)
#'
#' @export
#'
subsetOverlapProbMNk <- function(intMN, intAN, intBM, k){
    exponents <- vSum(vSubsetOverlapProbMNkNumerator(intMN, intAN, intBM, k),
                      -1 * vChoose(intMN, intBM))
    primes <- generate_n_primes(length(exponents))
    return(powerProduct(primes, exponents))
}

#' Compute the probability that two subsets of sets M and N intersect in
#' at least k points
#'
#' This function computes the probability that two subsets of sets M and N
#' intersect in at least k points.
#'
#' @inheritParams vSubsetOverlapProbMNkNumerator
#'
#' @return The probability that two subsets of sets M and N intersect in
#' at least k points.
#'
#' @examples
#' subsetOverlapPvalMNk(300, 23, 24, 6)
#'
#' @export
#'
subsetOverlapPvalMNk <- function(intMN, intAN, intBM, k){
    denom <- -1 * vChoose(intMN, intBM)
    exponents <- do.call(vSum,
                         lapply(seq(k, min(intAN, intBM)),
                                function(i) vSum(
                                    vSubsetOverlapProbMNkNumerator(
                                        intMN, intAN, intBM, i), denom)))
    primes <- generate_n_primes(length(exponents))
    return(powerProduct(primes, exponents))
}
