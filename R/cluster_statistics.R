#' @export

CountCCRSA <- function(markerNames, geneSets){
  clustersCCRSA <- lapply(markerNames, function(x) intersect(x, geneSets$Lit[[14]]))
  exclusiveCCRSA <- lapply(1:length(markerNames), function(x)
    setdiff(clustersCCRSA[[x]], Reduce(union, clustersCCRSA[setdiff(1:length(markerNames), x)])))
  return(unlist(mapply(function(x, y) c(length(x) - length(y), length(y)), clustersCCRSA, exclusiveCCRSA,
                       SIMPLIFY = F)))

}

ExclusiveClusterMarkers <- function(seuratObj, markers, message = T){
  counts <- sapply(1 + as.integer(levels(seuratObj)), function(x) length(setdiff(rownames(markers[[x]]), Reduce(union, lapply(markers[-x], rownames)))))
  if (message)return(MessageVector(counts))
  return(counts)
}

ExclusiveClusterMarkersPerc <- function(seuratObj, markers, message = T){
  counts <- ExclusiveClusterMarkers(seuratObj, markers, F)
  percentages <- round(counts/sum(counts) * 100, 2)
  if (message)return(MessageVector(str_c(percentages, "%")))
  return(percentages)
}

ClusterAvgExp <- function(seuratObj, features = rownames(seuratObj)){
  df <- data.frame(avgExp = AverageExpression(seuratObj, assay = "SCT", features = features, group.by = "seurat_clusters"))
  rownames(df) <- features
  colnames(df) <- levels(unique(seuratObj@meta.data[["seurat_clusters"]]))
  df$rowmax <- sapply(features, function(x) as.integer(colnames(df)[which.max(df[x,])]))
  df$rowmin <- sapply(features, function(x) as.integer(colnames(df)[which.min(df[x,-length(colnames(df))])]))
  return(df)
}

ClusterExtremaExpression <- function(seuratAvgExp, column, message = T){
  counts <- dplyr::count(seuratAvgExp, {{column}})$n
  if (message)return(MessageVector(counts))
  return(counts)
}

ClusterMaxExpression <- function(seuratAvgExp, message = T)
  return(ClusterExtremaExpression(seuratAvgExp, rowmax, message))

ClusterMinExpression <- function(seuratAvgExp, message = T)
  return(ClusterExtremaExpression(seuratAvgExp, rowmin, message))

ClusterExtremaExpressionPerc <- function(seuratAvgExp, column, message = T){
  percentages <- round(dplyr::count(seuratAvgExp, {{column}})$n/length(rownames(seuratAvgExp)) * 100, 2)
  if (message)return(MessageVector(str_c(percentages, "%")))
  return(percentages)
}

ClusterMaxExpressionPerc <- function(seuratAvgExp, message = T)
  return(ClusterExtremaExpressionPerc(seuratAvgExp, rowmax, message))

ClusterMinExpressionPerc <- function(seuratAvgExp, message = T)
  return(ClusterExtremaExpressionPerc(seuratAvgExp, rowmin, message))
