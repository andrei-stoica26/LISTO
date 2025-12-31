#' @importFrom stats dhyper median p.adjust phyper
#'
NULL


#' Perform multiple testing correction on a vector of p-values
#'
#' This function performs multiple testing correction on a vector of p-values.
#'
#' @param pvals A numeric vector.
#' @param mtMethod Multiple testing correction method. Choices are
#' 'BY' (default) 'holm', hochberg', hommel', 'bonferroni', 'BH',  'fdr' and
#' 'none'.
#' @param mtStat A statistics to be optionally computed. Choices are 'identity'
#' (no statistics will be computed and the adjusted p-values will be returned
#' as such), 'median', 'mean', 'max' and 'min'.
#' @param nComp Number of comparisons. In most situations, this parameter
#' should not be changed.
#'
#' @return A numeric vector.
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

    mtMethod <- match.arg(mtMethod, c('BY', 'holm', 'hochberg', 'hommel',
                                      'bonferroni', 'BH', 'fdr', 'none'))
    mtStat <- match.arg(mtStat, c('identity', 'median', 'mean', 'max', 'min'))
    statFun <- eval(as.name(mtStat))
    return(statFun(p.adjust(pvals, mtMethod, nComp)))
}

#' Perform multiple testing correction on a data frame
#'
#' This function orders a data frame based on a column of p-values, performs
#' multiple testing correction on the column, and filters the data-frame
#' based on the adjusted p-values.
#'
#' @inheritParams mtCorrectV
#' @param df A data frame with a p-values columnn.
#' @param colStr Name of the column of p-values.
#' @param newColStr Name of the column of adjusted p-values that will be
#' created.
#' @param doOrder Whether to increasingly order the data frame based on the
#' adjusted p-values.
#' @param doFilter Whether to filter the data frame based on the adjusted
#' p-values.
#' @param pvalThr p-value threshold used for filtering. Ignored
#' if \code{doFilter} is \code{FALSE}.
#'
#' @param ... Additional arguments passed to the multiple testing correction
#' method.
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
                        doOrder = TRUE,
                        doFilter = TRUE,
                        pvalThr = 0.05,
                        ...){
    mtMethod <- match.arg(mtMethod, c('BY', 'holm', 'hochberg', 'hommel',
                                      'bonferroni', 'BH', 'fdr', 'none'))
    df[[newColStr]] <- p.adjust(df[[colStr]], mtMethod, ...)
    if (doOrder)
        df <- df[order(df[[newColStr]]), ]
    if (doFilter)
        df <- df[df[, newColStr] < pvalThr, ]
    return(df)
}
