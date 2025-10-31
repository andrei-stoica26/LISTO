#' @importFrom stats setNames
#'
NULL

#' Extract two marker data frames from two named marker lists
#'
#' This function receives two named marker lists as input (they can be
#' identical) and extracts a specified marker data frame from each. The two
#' extracted marker data frames are grouped in a named list.
#'
#' @param markerList1 A named list of marker data frames.
#' @param markerList2 A named list of marker data frames.
#' @param name1 Name of the marker data frame that will be selected from
#'  \code{markerList1}.
#' @param name2 Name of the marker data frame that will be selected from
#'  \code{markerList2}.
#'
#' @return A named list of two marker data frames.
#'
#' @export
#'
markerDFPair <- function(markerList1, markerList2, name1, name2)
    return(setNames(list(markerList1[[name1]], markerList2[[name2]]),
                    c(name1, name2)))

#' Join two marker data frames
#'
#' This function joins two marker data frames.
#'
#' @param df1 A marker data frame.
#' @param df2 A marker data frame.
#' @param joinCol Column based on which the marker data frames will be joined.
#'
#' @return A data frame with the shared markers as rows and the join column
#' from each of the two marker data frames as columns.
#'
#' @export
#'
sharedMarkers <- function(df1, df2, joinCol = 'avg_log2FC'){
    shared <- intersect(rownames(df1), rownames(df2))
    df <- data.frame(cbind(df1[shared, joinCol],
                           df2[shared, joinCol]))
    rownames(df) <- shared
    colnames(df) <- paste0(joinCol, '_', c(1, 2))
    return(df)
}
