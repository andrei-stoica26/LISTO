#' Generate cutoffs for assessing marker overlaps
#'
#' This function generates cutoffs for assessing marker overlaps.
#'
#' @param obj1 A data frame with at least one numeric column.
#' @param obj2 A data frame with at least one numeric column.
#' @param col Numeric column used for ranking items.
#' @param isHighTop Whether higher values in the numeric column correspond to
#' top-ranked items.
#' @param extraCutoff Cutoff to be placed at one end of the cutoff list.
#' If \code{isHighTop} is \code{TRUE}, it must be lower than all the cutoffs
#' present in the numeric column. If \code{isHighTop} is \code{FALSE}, it must
#' be higher that all cutoffs present in the numeric column.
#' @param maxNCutoffs Maximum number of cutoffs. If the two input data frames
#' contain more cutoffs than this value, only \code{maxNCutoffs} linearly
#' spaced cutoffs will be selected from the original cutoff list.
#' @param verbose Whether the output should be verbose.
#'
#' @return A numeric vector.
#'
#' @keywords internal
#'
generateCutoffs <- function(obj1,
                            obj2,
                            col,
                            isHighTop = TRUE,
                            extraCutoff = 0,
                            maxNCutoffs = 500,
                            verbose = FALSE){

    values1 <- obj1[[col]]
    values2 <- obj2[[col]]
    cutoffs <- unique(c(values1, values2))
    if (isHighTop)
        cutoffs <- cutoffs[cutoffs < min(max(values1), max(values2))] else
            cutoffs <- cutoffs[cutoffs > max(min(values1), min(values2))]
    cutoffs <- c(cutoffs, extraCutoff)
    cutoffs <- sort(cutoffs, decreasing=isHighTop)
    nCutoffs <- length(cutoffs)
    if (nCutoffs > maxNCutoffs){
        if(verbose)
            message('Too many cutoffs found in the input data frames. Only ',
                    maxNCutoffs, ' will be used')
        cutoffs <- cutoffs[seq(1, nCutoffs, length.out=maxNCutoffs)]
    }
    return(cutoffs)
}
