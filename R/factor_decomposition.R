#' Generate the prime factor decomposition of n factorial.
#'
#' This function generates the prime factor decomposition of n factorial.
#'
#' @param n A positive integer.
#'
#' @return A vector in which positions represent prime numbers (that is, the
#' first position corresponds to 2, the second position corresponds to 3,
#' the third position corresponds to 5, etc.) and values
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

#' Compute the prime factor decomposition of the binomial coefficient
#'
#' This function computes the prime factor decomposition of the
#' binomial coefficient.
#'
#' @param n Total number of elements.
#' @param k Number of selected elements.
#'
#' @return The product of the primes raised to the exponents.
#'
#' @examples
#' powerProduct(c(2, 3, 5), c(4, 2, 6))
#'
#' @noRd
#'
powerProduct <- function(primes, exponents)
    return(prod(mapply(function(x, y) x ^ y, primes, exponents)))
