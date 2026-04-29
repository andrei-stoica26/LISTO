#' @importFrom stats dhyper median p.adjust phyper
#'
NULL

#' Perform Bonferroni correction on a vector of p-values
#'
#' This function performs Bonferroni correction on a vector of p-values.
#'
#' @details This function is implemented in order to allow users to perform
#' Bonferroni correction with fewer comparisons than the number of elements
#' in the vector, which is normally disallowed by the \code{p.adjust} function
#' from \code{stats}. A use case is correcting the Seurat markers returned
#' for each identity class by setting the number of comparisons to the number
#' of classes.
#'
#' @param pvals A numeric vector.
#' @param nComp Number of comparisons. In most situations, this parameter
#' should not be changed.
#'
#' @return Bonferroni-corrected p-values
#'
#' @keywords internal
#'
bfCorrectV <- function(pvals, nComp)
    return(vapply(pvals, function(x) min (x * nComp, 1), numeric(1)))


#' Helper function for multiple comparison testing
#'
#' This function is a helper for multiple comparison testing.
#'
#' @inheritParams bfCorrectV
#' @param mtMethod Multiple testing correction method. Choices are
#' 'BY' (default), 'holm', hochberg', hommel', 'bonferroni', 'BH',  'fdr' and
#' 'none'.
#'
#' @return Adjusted p-values.
#'
#' @keywords internal
#'
mtCorrectHelper <- function(pvals, mtMethod, nComp){
    if (nComp >= length(pvals))
        return(p.adjust(pvals, mtMethod, nComp))

    if (mtMethod == 'bonferroni')
        return(bfCorrectV(pvals, nComp))

    stop('`nComp` can be greater than the number of rows ',
         'only if `mtMethod` is set to bonferroni.')
}

#' Perform multiple testing correction on a vector of p-values
#'
#' This function performs multiple testing correction on a vector of p-values.
#'
#' @inheritParams mtCorrectHelper
#' @param mtStat A statistics to be optionally computed. Choices are 'identity'
#' (no statistics will be computed and the adjusted p-values will be returned
#' as such), 'median', 'mean', 'max' and 'min'.
#'
#' @return If \code{mtStat} is 'identity' (as default), a numeric vector of
#' p-values corrected for multiple testing. Otherwise, a statistic based on
#' these corrected p-values defined by \code{mtStat}.
#'
#' @examples
#' pvals <- c(0.032, 0.001, 0.0045, 0.051, 0.048)
#' mtCorrectV(pvals)
#'
#' @export
#'
mtCorrectV <- function(pvals,
                       mtMethod = c('BY', 'holm', 'hochberg', 'hommel',
                                    'bonferroni', 'BH', 'fdr', 'none'),
                       mtStat = c('identity', 'median', 'mean', 'max', 'min'),
                       nComp = length(pvals)){

    mtMethod <- match.arg(mtMethod)
    mtStat <- match.arg(mtStat)
    statFun <- eval(as.name(mtStat))
    return(statFun(mtCorrectHelper(pvals, mtMethod, nComp)))
}

#' Perform multiple testing correction on a data frame column
#'
#' This function orders a data frame based on a column of p-values, performs
#' multiple testing on the column, and filters the data-frame based on it.
#'
#' @inheritParams mtCorrectV
#' @param df A data frame with a p-values column.
#' @param colStr Name of the column of p-values.
#' @param newColStr Name of the column of adjusted p-values that will be
#' created.
#' @param doOrder Whether to increasingly order the data frame based on the
#' adjusted p-values.
#' @param pvalThr p-value threshold used for filtering. If \code{NULL}, no
#' filtering will be performed.
#'
#' @return A data frame with the p-value column corrected for multiple testing.
#'
#' @examples
#' df <- data.frame(elem = c('A', 'B', 'C', 'D', 'E'),
#' pval = c(0.032, 0.001, 0.0045, 0.051, 0.048))
#' mtCorrectDF(df)
#'
#' @export
#'
mtCorrectDF <- function(df,
                        mtMethod = c('BY', 'holm', 'hochberg', 'hommel',
                                     'bonferroni', 'BH', 'fdr', 'none'),
                        colStr = 'pval',
                        newColStr = 'pvalAdj',
                        pvalThr = 0.05,
                        doOrder = TRUE,
                        nComp = nrow(df)){
    mtMethod <- match.arg(mtMethod)
    df[[newColStr]] <- mtCorrectHelper(df[[colStr]], mtMethod, nComp)
    if (!is.null(pvalThr))
        df <- df[df[, newColStr] < pvalThr, ]
    if (doOrder)
        df <- df[order(df[[newColStr]]), ]
    return(df)
}

