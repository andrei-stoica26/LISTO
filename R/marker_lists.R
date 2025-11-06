#' @importFrom Seurat FindMarkers
#'
NULL

#' Find markers for Seurat identity class
#'
#' This function finds upregulated or downregulated markers for Seurat object
#' for a given identity class and performs an additional Bonferroni correction
#' for multiple testing.
#'
#' @inheritParams allGroups
#' @inheritParams allComplements
#' @param invert Whether to compute downregulated markers.
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
#' @export
#'
buildMarkerList <- function(seuratObj,
                            idClass = 'seurat_clusters',
                            invert = FALSE,
                            logFCThr = 0,
                            minPct = 0,
                            minPctRatio = 0,
                            ids1 = allGroups(seuratObj, idClass),
                            ids2 = allComplements(seuratObj, idClass, ids1),
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
                                   nTests=length(originalIds1))
            if(nrow(markers)){
                markers$pct.ratio <- markers$pct.1 / markers$pct.2
                markers <- markers[markers$pct.ratio >= minPctRatio, ]
            }

        }
        gc()
        return(markers)
    }, labels1, labels2, ids1, ids2, SIMPLIFY=FALSE)
    names(markerList) <- originalIds1
    markerList <- markerList[vapply(markerList, function(x) nrow(x) > 0,
                                    logical(1))]
    return(markerList)
}




