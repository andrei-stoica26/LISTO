#' @importFrom Seurat FindMarkers
#'
NULL

#' Find markers for a Seurat identity group
#'
#' This function finds markers for a Seurat identity group and performs an
#' additional Bonferroni correction for multiple testing.
#'
#' @inheritParams allComplements
#' @param id1 Selected group.
#' @param id2 Group selected for comparison.
#' @param logFCThr Fold change threshold for testing.
#' @param minPct The minimum fraction of in-cluster cells in which tested
#' genes need to be expressed.
#' @param minPctRatio The minimum ratio of in-cluster cells over out-cluster
#' cells in which a retained gene must be expressed.
#' @param ids2 Selected class groups used for comparison. Ignored
#' if \code{invert} is \code{TRUE}.
#' @param ... Additional arguments passed to \code{Seurat::FindMarkers}.
#'
#' @return A list of marker data frames.
#'
#' @keywords internal
#'
buildMarkerListCore <- function(seuratObj,
                                idClass = 'seurat_clusters',
                                id1 = NULL,
                                id2 = NULL,
                                logFCThr = 0,
                                minPct = 0,
                                minPctRatio = 0,
                                nTests = 1,
                                ...){

    markers <- FindMarkers(seuratObj,
                           group.by=idClass,
                           ident.1=id1,
                           ident.2=id2,
                           only.pos=TRUE,
                           logfc.threshold=logFCThr,
                           min.pct=minPct,
                           ...)
    if (nrow(markers)){
        markers <- mtCorrectDF(markers, 'bf', colStr='p_val_adj',
                                   nTests=nTests)
        if(nrow(markers)){
            markers$pct.ratio <- markers$pct.1 / markers$pct.2
            markers <- markers[markers$pct.ratio >= minPctRatio, ]
            }
    }
    gc()
    return(markers)
}

#' Find markers for Seurat identity groups
#'
#' This function finds upregulated or downregulated markers for Seurat object
#' for groups belonging to a given identity class and performs an additional
#' Bonferroni correction for multiple testing.
#'
#' @inheritParams allComplements
#' @inheritParams buildMarkerListCore
#' @param ids2 Selected class groups used for comparison. Ignored
#' if \code{invert} is \code{TRUE}.
#' @param invert Whether to compute downregulated markers rather than
#' upregulated ones.
#'
#' @return A list of marker data frames.
#'
#' @export
#'
buildMarkerList <- function(seuratObj,
                            idClass = 'seurat_clusters',
                            ids1 = allGroups(seuratObj, idClass),
                            ids2 = allComplements(seuratObj, idClass, ids1),
                            invert = FALSE,
                            logFCThr = 0,
                            minPct = 0,
                            minPctRatio = 0,
                            ...){
    originalIds1 <- ids1
    originalIds2 <- ids2
    if (invert){
        ids1 <- ids2
        ids2 <- originalIds1
        markerType <- 'downregulated'
    } else
        markerType <- 'upregulated'

    labels1 <- listToChar(originalIds1)
    labels2 <- listToChar(originalIds2)

    markerList <- mapply(function(label1, label2, id1, id2) {
        message('Finding ', markerType, ' markers for ', label1, '
                vs. ', label2, ' (', idClass, ')...')
        markers <- buildMarkerListCore(seuratObj, idClass, id1, id2,
                                       logFCThr, minPct, minPctRatio,
                                       length(labels1), ...)
        return(markers)
        }, labels1, labels2, ids1, ids2, SIMPLIFY=FALSE)
    names(markerList) <- originalIds1
    markerList <- markerList[vapply(markerList, function(x) nrow(x) > 0,
                                    logical(1))]
    return(markerList)
}

#' Find pairs of upregulated and downregulated markers for Seurat identity class
#'
#' This function finds pairs of upregulated or downregulated markers for
#' input groups belongign to a given identity class and performs an additional
#' Bonferroni correction for multiple testing.
#'
#' @inheritParams buildMarkerList
#'
#' @return A list of pairs of marker data frames.
#'
#' @export
#'
buildPairedMarkerList <- function(seuratObj,
                                  idClass = 'seurat_clusters',
                                  logFCThr = 0,
                                  minPct = 0,
                                  minPctRatio = 0,
                                  ids1 = allGroups(seuratObj, idClass),
                                  ids2 = allComplements(seuratObj, idClass,
                                                        ids1),
                                  ...){
    originalIds1 <- ids1
    originalIds2 <- ids2

    labels1 <- listToChar(ids1)
    labels2 <- listToChar(ids2)

    markerList <- mapply(function(label1, label2, id1, id2) {
        message('Finding upregulated markers for ', label1, '
                vs. ', label2, ' (', idClass, ')...')
        markers1 <- buildMarkerListCore(seuratObj, idClass, id1, id2,
                                        logFCThr, minPct, minPctRatio,
                                        length(labels1), ...)
        message('Finding upregulated markers for ', label2, '
                vs. ', label1, ' (', idClass, ')...')
        markers2 <- buildMarkerListCore(seuratObj, idClass, id2, id1,
                                        logFCThr, minPct, minPctRatio,
                                        length(labels1), ...)
        markerPairNames <- c(paste0(label1, ' vs. ', label2),
                             paste0(label2, ' vs. ', label1))
        return(setNames(list(markers1, markers2), markerPairNames))
    }, labels1, labels2, ids1, ids2, SIMPLIFY=FALSE)
    names(markerList) <- paste0(labels1, '_', labels2)
    return(markerList)
}
