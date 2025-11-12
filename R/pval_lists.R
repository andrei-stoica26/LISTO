#' Assess the overlap of two lists of marker data frames.
#'
#' This function assesses the overlap of two lists of marker data frames.
#'
#' @inheritParams pvalObjects
#' @param list1 List of data frames.
#' @param list2 List of data frames.
#' @param mtMethod Multiple testing correction method.
#' @param ... Additional arguments passed to \code{mtCorrectDF}.
#'
#' @return A data frame.
#'
#' @export
#'
pvalLists <- function(list1,
                      list2,
                      col,
                      allItems1,
                      allItems2 = NULL,
                      nCores = 1,
                      isHighTop = TRUE,
                      extraCutoff = 0,
                      maxNCutoffs = 500,
                      mtMethod = c('BY', 'holm', 'hochberg',
                                   'hommel', 'bonferroni', 'BH',
                                   'fdr', 'none'),
                      verbose = FALSE,
                      ...){

    mtMethod <- match.arg(mtMethod, c('BY', 'holm', 'hochberg',
                                      'hommel', 'bonferroni', 'BH',
                                      'fdr', 'none'))

    df <- expand.grid(names(list1), names(list2))
    colnames(df) <- c('Group1', 'Group2')

    if (is.null(allItems2) & nCores > 1)
        message('Parallelization is not supported for single-dataset ',
                'overlap assessments. `nCores` will be ignored.')

    df$pval <- vapply(seq(nrow(df)),
                      function(i){
                          names1 <- df[i, 1]
                          names2 <- df[i, 2]
                          message('Assessing overlap between sets: ',
                                names1, ' and ', names2, '...')
                          pvalObjects(
                              list1[[names1]],
                              list2[[names2]],
                              col,
                              allItems1,
                              allItems2,
                              nCores,
                              isHighTop,
                              extraCutoff,
                              maxNCutoffs,
                              mtMethod,
                              verbose)
                          }, numeric(1))

    df <- mtCorrectDF(df, mtMethod, ...)
    return(df)
}
