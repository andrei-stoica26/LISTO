#' @importFrom henna riverPlot
#'
NULL

#' Plot representation data frame
#'
#' This function plots representation data frame as an alluvial plot.
#'
#' @inheritParams prepAlluvial
#' @param ... Additional parameters passed to \code{henna::riverPlot}
#'
#' @return A ggplot object
#'
#' @export
#'
pvalRiverPlot <- function(df,
                          pvalCol = 'pvalAdj',
                          colIndices = c(1, 2),
                          weightExp = 1/2,
                          pvalOffset = 1e-317,
                          ...){
    resDF <- prepAlluvial(df, pvalCol, colIndices, weightExp, pvalOffset)
    p <- riverPlot(resDF, ...)
    return(p)
}
