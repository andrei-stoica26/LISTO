#' @importFrom henna hullPlot
#'
NULL

#' Find markers for Seurat identity class
#'
#' This function finds upregulated or downregulate markers for Seurat object
#' for a given identity class and performs an additional Bonferroni correction for
#' multiple testing.
#'
#' @inheritParams markerListPhyper
#' @param name1 Name of first marker data frame.
#' @param name2 Name of second marker data frame.
#' @inheritParams sharedMarkers
#' @param thresh1 Join column threshold for the first marker data frame.
#' @param thresh2 Join column threshold for the second marker data frame.
#' @param title Plot title
#' @param ... Additional arguments passed to \code{henna::hullPlot}.
#'
#' @return A list of marker data frames.
#'
#' @export
#'
sharedMarkersPlot <- function(markerList1,
                              markerList2,
                              name1,
                              name2,
                              joinColumn = 'avg_log2FC',
                              thresh1 = 1.5,
                              thresh2 = 1.5,
                              title = paste0('Shared markers - ', name1,
                                        ' and ', name2),
                              ...){
    sharedDF <- sharedMarkers(markerList1[[name1]],
                              markerList2[[name2]])
    labelDF <- sharedDF[sharedDF[, 1] > thresh1 & sharedDF[, 2] > thresh2, ]
    p <- hullPlot(sharedDF,
                  title,
                  xInt=thresh1,
                  yInt=thresh2,
                  xLab=paste0(joinColumn, ' (', name1, ')'),
                  yLab=paste0(joinColumn, ' (', name2, ')'),
                  legendLabs=as.factor(c('Non-top markers',
                               'Shared top markers',
                               paste0('Top markers only for ', name2),
                               paste0('Top markers only for ', name1))),
                  labelDF = labelDF,
                  ...
                  )
    return(p)
}
