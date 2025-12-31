#' Build a Seurat marker list ready to be used by LISTO
#'
#' This function builds a Seurat marker list ready to be used by LISTO.
#' Requires Seurat (not automatically installed with LISTO).
#'
#' @param seuratObj A Seurat object.
#' @param col Seurat metadata column used for grouping.
#' @param logFCThr Fold change threshold for testing.
#' @param minPct The minimum fraction of in-cluster cells in which tested
#' genes need to be expressed.
#' @param ... Additional arguments passed to \code{Seurat::FindMarkers}.
#'
#' @export
#'
buildSeuratMarkerList <- function(seuratObj,
                                  col,
                                  logFCThr = 1,
                                  minPct = 0.2,
                                  ...){
    if (requireNamespace('Seurat', quietly=TRUE)){
        groups <- unique(seuratObj[[]][[col]])
        markers <- lapply(groups,
                          function(x){
                              message('Computing markers for ', x, '...')
                              Seurat::FindMarkers(seuratObj,
                                                  group.by=col,
                                                  ident.1=x,
                                                  only.pos=TRUE,
                                                  min.pct=minPct,
                                                  logfc.threshold=logFCThr,
                                                  ...)
                               })
        names(markers) <- groups
        return(markers)
    } else
        print('Install Seurat to use `buildSeuratMarkerList`.')
}
