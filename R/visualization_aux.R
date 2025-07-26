#' Join two marker data frames
#'
#' This function joins two marker data frames.
#'
#' @param markerSet1 A marker data frame.
#' @param markerSet2 A marker data frame.
#' @param joinCol Join column.
#'
#' @return A data frame with the shared markers as rows and the join column
#' from each of the two marker data frames as columns.
#'
#' @export
#'
sharedMarkers <- function(markerSet1, markerSet2, joinCol = 'avg_log2FC'){
    shared <- intersect(rownames(markerSet1), rownames(markerSet2))
    df <- data.frame(cbind(markerSet1[shared, joinCol],
                           markerSet2[shared, joinCol]))
    rownames(df) <- shared
    colnames(df) <- paste0(joinCol, '_', c(1, 2))
    return(df)
}
