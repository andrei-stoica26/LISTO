#' Check if numCol is valid input for an object.
#'
#' This function checks if \code{numCol} is valid input for an object.
#'
#' @inheritParams getObjectValues
#'
#' @return Nothing. This function is called for its side effect.
#'
#' @noRd
#'
checkNumCol <- function(obj, numCol){
    if(is(obj, 'data.frame')){
        if (!numCol %in% colnames(obj))
            stop('All data frames must contain a `numCol` column.')
        if (!is(obj[, numCol], 'numeric'))
            stop('The `numCol` column must be numeric in all data frames.')
    }

}

#' Check if numCol is valid input for a list of objects
#'
#' This function checks if \code{numCol} is valid input for a list of objects.
#'
#' @inheritParams getObjectValues
#' @param objs A list containing data frame with a numeric column
#' or character vectors.
#'
#' @return Nothing. This function is called for its side effect.
#'
#' @noRd
#'
checkNumColAll <- function(objs, numCol)
    for (obj in objs)
        checkNumCol(obj, numCol)
