getObjectValues <- function(obj, numCol = NULL, isHighTop = TRUE){
    if (is.null(numCol) | is(obj, 'character')){
        if(isHighTop)
            return (-Inf)
        return(Inf)
    }
    return(obj[[numCol]])
}

#' @param isHighTop Whether higher values in the numeric column correspond to
#' top-ranked items.
#' @param extraCutoff Cutoff to be placed at one end of the cutoff list.
#' If \code{isHighTop} is \code{TRUE}, it must be lower than all the cutoffs
#' present in the numeric column. If \code{isHighTop} is \code{FALSE}, it must
#' be higher that all cutoffs present in the numeric column.
#' @param maxCutoffs Maximum number of cutoffs. If the input data frames
#' contain more cutoffs than this value, only \code{maxCutoffs} linearly
#' spaced cutoffs will be selected from the original cutoff list.
#'
#' @return A numeric vector.
#'
#' @keywords internal
#'
generateCutoffs <- function(obj1,
                            obj2,
                            obj3 = NULL,
                            numCol = NULL,
                            isHighTop = TRUE,
                            maxCutoffs = 500){

    values1 <- getObjectValues(obj1, numCol, isHighTop)
    values2 <- getObjectValues(obj2, numCol, isHighTop)
    if(is.null(obj3)){
        values3 <- values2
        cutoffs <- unique(c(values1, values2))
    } else {
        values3 <- getObjectValues(obj3, numCol, isHighTop)
        cutoffs <- unique(c(values1, values2, values3))
    }

    if (isHighTop){
        bound <- min(max(values1), max(values2), max(values3))
        cutoffs <- cutoffs[cutoffs < bound]
    } else{
        bound <- max(min(values1), min(values2), min(values3))
        cutoffs <- cutoffs[cutoffs > bound]
    }

    cutoffs <- sort(cutoffs, decreasing=isHighTop)
    extraCutoff <- (1 - 2 * as.integer(isHighTop)) * Inf
    cutoffs <- unique(c(cutoffs, extraCutoff))
    nCutoffs <- length(cutoffs)
    if (nCutoffs > maxCutoffs){
        message(nCutoff, 'cutoffs found in the input data frames. Only ',
                maxCutoffs, ' will be used. To change this behavior, set a ',
                'higher value to `maxCutoffs`.')
        cutoffs <- cutoffs[seq(1, nCutoffs, length.out=maxCutoffs)]
    }
    return(cutoffs)
}
