#' @importFrom Seurat FindMarkers
#'
NULL

#' Find upregulated markers for Seurat identity class
#'
#' This function finds upregulated markers for Seurat object for a given
#' identity class and performs an additional Bonferroni correction for
#' multiple testing.
#'
#' @param seuratObj A Seurat object.
#' @param idClass Identity class.
#' @param logfcThreshold Fold change threshold for testing.
#' @param identities Selected class groups.
#'
#' @return A list of marker data frames.
#'
#' @export
#'
findUpMarkers <- function(seuratObj,
                          idClass = 'seurat_clusters',
                          logfcThreshold = 0,
                          identities = sort(unique(seuratObj[[idClass]]))){
    markerList <- lapply(identities, function(x){
        message('Finding upregulated markers for identity class ', x, '...')
        markers <- FindMarkers(seuratObj,
                               group.by=idClass,
                               ident.1=x,
                               only.pos=TRUE,
                               logfc.threshold=logfcThreshold,
                               min.pct=0,
                               densify=TRUE)
        return(hammers::bfCorrectDF(markers, length(identities)))
    })
    names(markerList) <- identities
    return(markerList)
}

#' Find downregulated markers for Seurat identity class
#'
#' This function finds downregulated markers for Seurat object for a given
#' identity class and performs an additional Bonferroni correction for
#' multiple testing.
#'
#' @inheritParams findUpMarkers
#'
#' @export
#'
findDownMarkers <- function(seuratObj,
                          idClass = 'seurat_clusters',
                          logfcThreshold = 0,
                          identities = sort(unique(seuratObj[[]][[idClass]]))){
    allIdentities <- sort(unique(seuratObj[[]][[idClass]]))
    markerList <- lapply(identities, function(x) {
        message('Finding downregulated markers for identity class ', x, '...')
        markers <- FindMarkers(seuratObj,
                               group.by=idClass,
                               ident.1=setdiff(allIdentities, x),
                               only.pos=TRUE,
                               logfc.threshold=logfcThreshold,
                               min.pct=0,
                               densify=TRUE)
        return(hammers::bfCorrectDF(markers, length(identities)))
    })
    names(markerList) <- identities
    return(markerList)
}
