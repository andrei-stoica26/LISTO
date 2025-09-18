#' @importFrom henna hullPlot
#'
NULL

#' Find markers for Seurat identity class
#'
#' This function finds upregulated or downregulate markers for Seurat object
#' for a given identity class and performs an additional Bonferroni correction
#' for multiple testing.
#'
#' @inheritParams markerListPhyper
#' @param name1 Name of first marker data frame.
#' @param name2 Name of second marker data frame.
#' @param markerNames Names of markers to be displayed on the plot.
#' If \code{NULL}, the markers will be chosen based on the thresholds.
#' @inheritParams sharedMarkers
#' @param thresh1 Join column threshold for the first marker data frame.
#' @param thresh2 Join column threshold for the second marker data frame.
#' @param title Plot title.
#' @param markersHullColor Color of the hull determined by \code{markers}.
#' Ignored if \code{markers} is \code{NULL}.
#' @param isNeg1 Whether the first marker set represents negative markers.
#' @param isNeg2 Whether the second marker set represents negative markers.
#' @param ... Additional arguments passed to \code{henna::hullPlot}.
#'
#' @return A list of marker data frames.
#'
#' @export
#'
sharedMarkersPlot <- function(markerList1,
                              name1,
                              name2,
                              markerList2 = NULL,
                              markerNames = NULL,
                              joinColumn = 'avg_log2FC',
                              thresh1 = 0,
                              thresh2 = 0,
                              title = 'Shared markers',
                              markersHullColor = 'purple',
                              isNeg1 = FALSE,
                              isNeg2 = FALSE,
                              ...){

    if(is.null(markerList2))
        markerList2 <- markerList1

    sharedDF <- sharedMarkers(markerList1[[name1]],
                              markerList2[[name2]])

    if(isNeg1)
        name1 <- paste0(name1, ' (downregulated)')
    if(isNeg2)
        name2 <- paste0(name2, ' (downregulated)')

    if (is.null(markerNames)){
        labelDF <- sharedDF[sharedDF[, 1] > thresh1 &
                                sharedDF[, 2] > thresh2, ]
        p <- hullPlot(sharedDF,
                      title,
                      xInt=thresh1,
                      yInt=thresh2,
                      xLab=paste0(joinColumn, ' (', name1, ')'),
                      yLab=paste0(joinColumn, ' (', name2, ')'),
                      legendLabs=as.factor(c('Non-top',
                                             'Shared',
                                             paste0('Top only for ', name2),
                                             paste0('Top only for ', name1))),
                      labelDF=labelDF,
                      ...
        )

    }else{
        foundMarkers <- intersect(rownames(sharedDF), markerNames)
        message(length(foundMarkers), ' among the ', length(markerNames),
                ' markers are shared by the ', name1, ' and ',
                name2, ' markers.')
        labelDF <- sharedDF[foundMarkers, ]
        p <- hullPlot(labelDF,
                      title,
                      xLab=paste0(joinColumn, ' (', name1, ')'),
                      yLab=paste0(joinColumn, ' (', name2, ')'),
                      legendLabs='Shared markers',
                      labelDF=labelDF,
                      palette=markersHullColor,
                      ...
        )
    }
    return(p)
}
