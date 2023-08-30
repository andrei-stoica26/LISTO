#Find markers of treatment condition selections
FindSelectionMarkers <- function(seuratObj, RownamesCML, firsts, seconds, nTests, logfcThreshold = 0, recorrectUMI = T)
{
  indices <- 1:length(RownamesCML)
  return(lapply(indices, function(x) {
    markers <- FindMarkers(seuratObj, group.by = "orig.ident", ident.1 = firsts[[x]],
                           ident.2 = seconds[[x]], only.pos = T, logfc.threshold = logfcThreshold,
                           min.pct = 0, densify = T, recorrect_umi = recorrectUMI)
    gc()
    return(BonferroniAndFilter(markers, nTests))
  }))
}
#' @export

#Find genes significantly upregulated in clusters
FindUpMarkers <- function(seuratObj, logfcThreshold = 0, indices = 0:(length(levels(seuratObj)) - 1))
  #Assumes clusters are ordered from 0 to maxNCluster
  return(lapply(indices, function(x) {
    markers <- FindMarkers(seuratObj, ident.1 = x, only.pos = T, logfc.threshold = logfcThreshold,
                           min.pct = 0, densify = T)
    gc()
    return(BonferroniAndFilter(markers, length(indices)))
  }))

#Find genes significantly downregulated in clusters
FindDownMarkers <- function(seuratObj, logfcThreshold = 0, indices = 0:(length(levels(seuratObj)) - 1))
  #Assumes clusters are ordered from 0 to maxNCluster
  return(lapply(indices, function(x) {
    markers <- FindMarkers(seuratObj, ident.1 = setdiff(indices, x), only.pos = T, logfc.threshold = logfcThreshold, min.pct = 0, densify = T)
    gc()
    return(BonferroniAndFilter(markers, length(indices)))
  }))

#Find markers of condition selections within clusters
FindICSMarkers <- function(seurats, RownamesCML, firsts, seconds, logfc.threshold = 0){
  return(invisible(lapply(1:length(seurats), function(i) {
    message(str_c("Calculating selection markers in Cluster ", i - 1))
    FindSelectionMarkers(seurats[[i]], RownamesCML, firsts, seconds, logfc.threshold,
                         nTests = length(RownamesCMLA13A) * length(seurats), recorrectUMI = F)
  })))
}

