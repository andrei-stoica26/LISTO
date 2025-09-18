#' @importFrom primes generate_primes
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
choose <- function(n, k)
    return(vSum(factorialPrimePowers(n),
                -factorialPrimePowers(k),
                -factorialPrimePowers(n - k)))
