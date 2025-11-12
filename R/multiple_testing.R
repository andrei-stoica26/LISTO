#' Perform multiple testing correction on a vector of p-values
#'
#' This function perform multiple testing correction on a vector of p-values.
#'
#' @param pvals A numeric vector.
#' @param mtMethod Multiple testing correction method. Choices are
#' Bonferroni ('bf'), Benjamini-Hochberg('bh'), and Benjamini-Yekutieli ('by').
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
                       mtMethod = c('holm', 'hochberg', 'hommel',
                                    'bonferroni', 'BH', 'BY',
                                    'fdr', 'none'),
                       mtStat = c('identity', 'median', 'mean', 'max', 'min'),
                       nComp = length(pvals)){

    mtMethod <- match.arg(mtMethod, c('holm', 'hochberg', 'hommel',
                                      'bonferroni', 'BH', 'BY',
                                      'fdr', 'none'))
    mtStat <- match.arg(mtStat, c('identity', 'median', 'mean', 'max', 'min'))
    statFun <- eval(as.name(mtStat))
    return(statFun(p.adjust(pvals, mtMethod, nComp)))
}

#' Perform multiple testing correction on a data frame column
#'
#' This function orders a data frame based on a column of p-values, performs
#' multiple testing on the column, and filters the data-frame based on it.
#'
#' @inheritParams mtCorrectV
#' @param df A data frame with a p-values columnn.
#' @param pvalThr p-value threshold.
#' @param colStr Name of the column of p-values.
#' @param newColStr Name of the column of adjusted p-values that will be
#' created.
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
                        mtMethod = c('holm', 'hochberg', 'hommel',
                                     'bonferroni', 'BH', 'BY',
                                     'fdr', 'none'),
                        pvalThr = 0.05,
                        colStr = 'pval',
                        newColStr = 'pvalAdj'){
    df <- df[order(df[[colStr]]), ]
    df[[newColStr]] <- mtCorrectV(df[[colStr]], mtMethod)
    df <- df[df[, newColStr] < pvalThr, ]
    return(df)
}
