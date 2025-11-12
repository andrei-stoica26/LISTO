#' @importFrom henna hullPlot volcanoPlot
#'
NULL

#' Plots marker shared between two marker data frames
#'
#' This function plots markers shared between two marker data frames.
#'
#' @details A wrapper around \code{henna::hullPlot}
#'
#' @inheritParams markerDFListOverlap
#' @inheritParams sharedMarkers
#' @param markerObj A list of two named marker data frames.
#' @param title Plot title.
#' @param labeledMarkers Names of markers to be displayed on the plot.
#' If \code{NULL}, the markers will be chosen based on the thresholds.
#' @param thresh1 Join column threshold for the first marker data frame.
#' @param thresh2 Join column threshold for the second marker data frame.
#' @param isNeg1 Whether the first marker set represents negative markers.
#' @param isNeg2 Whether the second marker set represents negative markers.
#' @param nameSuffix1 Suffix appended to \code{name1} on axes and legend.
#' @param nameSuffix2 Suffix appended to \code{name2} on axes and legend.
#' @param verbose Whether the output should specify how many of the markers in
#' \code{labeledMarkers} are top markers. Ignored if \code{labeledMarkers}
#' is NULL.
#' @param ... Additional arguments passed to \code{henna::hullPlot}.
#'
#' @return A ggplot object.
#'
#' @export
#'
sharedMarkersPlot <- function(markerObj,
                              title = NULL,
                              labeledMarkers = NULL,
                              joinCol = 'avg_log2FC',
                              thresh1 = NULL,
                              thresh2 = NULL,
                              isNeg1 = FALSE,
                              isNeg2 = FALSE,
                              nameSuffix1 = NULL,
                              nameSuffix2 = NULL,
                              verbose = TRUE,
                              ...){

    sharedDF <- sharedMarkers(markerObj[[1]], markerObj[[2]], joinCol)

    name1 <- paste0(names(markerObj)[1], nameSuffix1)
    name2 <- paste0(names(markerObj)[2], nameSuffix2)

    if(isNeg1)
        name1 <- paste0(name1, ' (downregulated)')
    if(isNeg2)
        name2 <- paste0(name2, ' (downregulated)')

    if (is.null(labeledMarkers) & verbose){
        topMarkers <- rownames(sharedDF[sharedDF[, 1] > thresh1 &
                                            sharedDF[, 2] > thresh2, ])
        foundMarkers <- intersect(topMarkers, labeledMarkers)
        message(length(foundMarkers), ' among the ', length(labeledMarkers),
                ' input markers are shared by the ', name1, ' and ',
                name2, ' markers.')
    }

    p <- hullPlot(sharedDF,
                  title,
                  xInt=thresh1,
                  yInt=thresh2,
                  xLab=paste0(joinCol, ' (', name1, ')'),
                  yLab=paste0(joinCol, ' (', name2, ')'),
                  legendLabs=as.factor(c('Non-top',
                                         'Shared',
                                         paste0('Top only for ', name2),
                                         paste0('Top only for ', name1))),
                  labeledPoints=labeledMarkers,
                  labXThr=thresh1,
                  labYThr=thresh2,
                  ...
                  )
    return(p)
}
