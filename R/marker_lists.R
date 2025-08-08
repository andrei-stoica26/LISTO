#' @importFrom Seurat FindMarkers
#'
NULL

#' Find markers for Seurat identity class
#'
#' This function finds upregulated or downregulate markers for Seurat object
#' for a given identity class and performs an additional Bonferroni correction for
#' multiple testing.
#'
#' @param seuratObj A Seurat object.
#' @param idClass Identity class.
#' @param invert Whether to compute downregulated markers.
#' @param logFCThr Fold change threshold for testing.
#' @param ids1 Selected class groups.
#' @param ids2 Selected class groups used for comparison. Ignored
#' if \code{invert} is \code{TRUE}.
#' @param ... Additional arguments passed to FindMarkers.
#'
#' @return A list of marker data frames.
#'
#' @export
#'
buildMarkerList <- function(seuratObj,
                            idClass = 'seurat_clusters',
                            invert = FALSE,
                            logFCThr = 0,
                            ids1 = sort(unique(seuratObj[[]][[idClass]])),
                            ids2 = NULL,
                            ...){
    originalIds1 <- ids1
    allIdentities <- sort(unique(seuratObj[[]][[idClass]]))
    diffs <- lapply(ids1, function(x) setdiff(allIdentities, x))
    if (invert){
        ids1 <- diffs
        ids2 <- originalIds1
        markerType <- 'downregulated'
    } else{
        if (is.null(ids2))
            ids2 <- diffs
        markerType <- 'upregulated'
    }
    markerList <- mapply(function(origId1, id1, id2) {
        message('Finding ', markerType, ' markers for identity class ',
                origId1, '...')
        markers <- FindMarkers(seuratObj,
                               group.by=idClass,
                               ident.1=id1,
                               ident.2=id2,
                               only.pos=TRUE,
                               logfc.threshold=logFCThr,
                               min.pct=0,
                               ...)
        if (nrow(markers)){
            markers <- bfCorrectDF(markers, length(originalIds1),
                                            colStr='p_val_adj')
            if(nrow(markers))
                markers$pct.ratio <- markers$pct.1 / markers$pct.2
        }
        gc()
        return(markers)
    }, originalIds1, ids1, ids2, SIMPLIFY=FALSE)
    names(markerList) <- originalIds1
    markerList <- markerList[vapply(markerList, function(x) nrow(x) > 0,
                                    logical(1))]
    return(markerList)
}
