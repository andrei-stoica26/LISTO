#' Assess the overlap of two lists of marker data frames.
#'
#' This function assesses the overlap of two lists of marker data frames.
#'
#' @inheritParams pvalObjects
#' @param list1 A list containing character vectors or data frames having one
#' numeric column.
#' @param list2 A list containing character vectors or data frames having one
#' numeric column.
#' @param list3 A list containing character vectors or data frames having one
#' numeric column.
#' @param universe1 The set from which the items corresponding to the elements
#' in \code{list1} are selected.
#' @param universe2 The set from which the items corresponding to the elements
#' in \code{list2} are selected.
#' @param mtMethod Multiple testing correction method.
#' @param ... Additional arguments passed to \code{mtCorrectDF}.
#'
#' @return A data frame.
#'
#' @export
#'
runLISTO <- function(list1,
                     list2,
                     list3 = NULL,
                     universe1,
                     universe2 = NULL,
                     numCol = NULL,
                     isHighTop = TRUE,
                     maxCutoffs = 500,
                     mtMethod = c('BY', 'holm', 'hochberg',
                                  'hommel', 'bonferroni', 'BH',
                                  'fdr', 'none'),
                     nCores = 1,
                     ...){
    mtMethod <- match.arg(mtMethod, c('BY', 'holm', 'hochberg',
                                      'hommel', 'bonferroni', 'BH',
                                      'fdr', 'none'))

    if (is.null(list3)){
        df <- expand.grid(names(list1), names(list2))
        if (is.null(universe2))
            type <- '2N' else
                type <- '2MN'
    } else {
        df <- expand.grid(names(list1), names(list2), names(list3))
        type <- '3N'
        if (!is.null(universe2))
            message('Three-way overlaps can be currently computed only for',
                    ' one universe. `universe2` will be ignored.')
    }
    colnames(df) <- paste0('Group', seq(ncol(df)))

    df$pval <- vapply(seq(nrow(df)),
                      function(i){
                          names1 <- df[i, 1]
                          names2 <- df[i, 2]
                          obj1 <- list1[[names1]]
                          obj2 <- list2[[names2]]
                          if(is.null(list3)){
                              obj3 <- NULL
                              message('Assessing overlap between sets ',
                                      names1, ' and ', names2, '...')
                          } else {
                              names3 <- df[i, 3]
                              obj3 <- list3[[names3]]
                              message('Assessing overlap between sets ',
                                      names1, ', ', names2, ' and ', names3,
                                      '...')
                          }
                          return(pvalObjects(obj1, obj2, obj3,
                                             universe1, universe2, numCol,
                                             isHighTop, maxCutoffs, mtMethod,
                                             nCores, type))
                      }, numeric(1))
    df <- mtCorrectDF(df, mtMethod, ...)
    return(df)
}
