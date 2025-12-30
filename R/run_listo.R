runLISTO <- function(list1,
                     list2,
                     list3 = NULL
                     universe1,
                     universe2 = NULL,
                     numCol = NULL,
                     isHighTop = TRUE,
                     maxCutoffs = 500,
                     mtMethod = c('BY', 'holm', 'hochberg',
                                  'hommel', 'bonferroni', 'BH',
                                  'fdr', 'none'),
                     nCores = 1){
    mtMethod <- match.arg(mtMethod, c('BY', 'holm', 'hochberg',
                                      'hommel', 'bonferroni', 'BH',
                                      'fdr', 'none'))
    if(is.null(list3))
        df <- expand.grid(names(list1), names(list2)) else
            df <- expand.grid(names(list1), names(list2), names(list3))
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
                          return(pvalObjects(obj1, obj2, obj3, universe1,
                                             universe2, numCol, isHighTop,
                                             maxCutoffs, mtMethod, nCores))
                      }, numeric(1))
    df <- mtCorrectDF(df, mtMethod, ...)
    return(df)
}
