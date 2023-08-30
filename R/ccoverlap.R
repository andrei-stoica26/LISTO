#' @export

#Overlaps between condition markers and cluster markers - NOT intracluster
CCOverlap <- function(seuratObj, clusterMarkers, selectionMarkers, RownamesCML, typeNames, nItems, isLog = T){
  clusters <- c()
  pvalues <- c()
  selections <- c()
  nGenes <- length(rownames(seuratObj))
  for (i in 1:length(levels(seuratObj$seurat_clusters)))
  {
    message(str_c("Calculating overlaps with cluster ", i - 1, "."))
    for (j in 1:length(RownamesCML)){
      message(str_c("Calculating overlaps with condition ", RownamesCML[j], "."))
      clusters <- c(clusters, typeNames[i])
      selections <- c(selections, RownamesCML[j])
      pvalues <- c(pvalues, TwoLogs(clusterMarkers[[i]], selectionMarkers[[j]], nItems, isLog))
    }
  }
  df <- data.frame(Cluster = clusters, Grouping = selections, pvalue = pvalues)
  return(BYCorrectDF(df))
}

CCPathwaysOverlap <- function(seuratObj, clusterDescriptions, selectionDescriptions, RownamesCML, typeNames, nItems)
  return(CCOverlap(seuratObj, clusterDescriptions, selectionDescriptions, RownamesCML, typeNames, nItems, isLog = F))

#Three way overlaps

ThreeWayCCOverlap <- function(thirdSet, clusterMarkers, selectionMarkers, RownamesCML, typeNames, N, ThreeFun){
  clusters <- c()
  pvalues <- c()
  selections <- c()
  for (i in 1:length(clusterMarkers))
  {
    message(str_c("Calculating overlaps with cluster ", i - 1, "."))
    for (j in 1:length(RownamesCML)){
      message(str_c("Calculating overlaps with condition ", RownamesCML[j], "."))
      clusters <- c(clusters, typeNames[i])
      selections <- c(selections, RownamesCML[j])
      pvalues <- c(pvalues, ThreeFun(thirdSet, clusterMarkers[[i]], selectionMarkers[[j]], N))
    }
  }
  df <- data.frame(Cluster = clusters, Grouping = selections, pvalue = pvalues)
  return(BYCorrectDF(df))
}

LitCCOverlap <- function(thirdSet, clusterMarkers, selectionMarkers, RownamesCML, typeNames, N)
  return(ThreeWayCCOverlap(thirdSet, clusterMarkers, selectionMarkers, RownamesCML, typeNames, N, OneLitTwoLogs))

LogCCOverlap <- function(thirdSet, clusterMarkers, selectionMarkers, RownamesCML, typeNames, N)
  return(ThreeWayCCOverlap(thirdSet, clusterMarkers, selectionMarkers, RownamesCML, typeNames, N, ThreeLogs))

AllLitCCOverlap <- function(allLit, clusterMarkers, selectionMarkers, RownamesCML, typeNames, N){
  return(lapply(allLit, function(x) LitCCOverlap(x, clusterMarkers, selectionMarkers, RownamesCML, typeNames, N)))
}

PrintAllLitCCOverlap <- function(CC, LitNames){
  invisible(lapply(1:length(LitNames), function(x){
    print(LitNames[x])
    PrintOverlapCC(CC[[x]])
  }))
}
