#' Prepare dataframe for alluvial plot
#'
#' This function extracts the relevant information from dataframe and adjusts
#' p-values to be used as weights for the alluvia.
#'
#' @param df A data frame.
#' @param pvalCol Name of p-value column to be used by the alluvial plot.
#' @param colIndices A vector respresenting the indices of the two categorical
#' columns from the data frame that will be used.
#' @param weightExp Exponent used in constructing weight from p-values.
#' @param pvalOffset Offset used to avoid zeros inside the logarithm function.
#'
#' @return A data frame with weight scores in lieu of p-values.
#'
#' @keywords internal
#'
prepAlluvial <- function(df,
                         pvalCol = 'pvalAdj',
                         colIndices = c(1, 2),
                         weightExp = 1/2,
                         pvalOffset = 1e-317){
    pvals <- sort(df[[pvalCol]])
    pvals[-1] <- -log(pvals[-1] + pvalOffset)
    if (pvals[1])
        pvals[1] <- -log(pvals[1]) else
            if (length(pvals) > 2)
                pvals[1] <- 2 * pvals[2] - pvals[3] else
                    pvals[1] <- 1
    resDF <- df[, colIndices]
    resDF$weight <- pvals ^ weightExp
    return(resDF)
}
