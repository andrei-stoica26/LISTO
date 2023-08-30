#' @export

ModifyER <- function(ER, nGroupings){
  ER@result$p.adjust <- ER@result$p.adjust * nGroupings
  ER@result <- subset(ER@result, p.adjust < 0.05)
  return(ER)
}

MarkersEnrich <- function(markerList, nGroupings, database)
  return(lapply(markerList, function(x){
    ER <- GenesER(rownames(x), database)
    return(ModifyER(ER, nGroupings))
  }))


ClusterGO <- function(seuratObj, markerList)
  return(MarkersEnrich(markerList, length(levels(seuratObj$seurat_clusters)), "enrichGO"))

ClusterKEGG <- function(seuratObj, markerList)
  return(MarkersEnrich(markerList, length(levels(seuratObj$seurat_clusters)), "enrichKEGG"))

ClusterWP<- function(seuratObj, markerList)
  return(MarkersEnrich(markerList, length(levels(seuratObj$seurat_clusters)), "enrichWP"))


SelectionGO <- function(RownamesCML, markerList)
  return(MarkersEnrich(markerList, length(RownamesCML), "enrichGO"))

SelectionKEGG <- function(RownamesCML, markerList)
  return(MarkersEnrich(markerList, length(RownamesCML), "enrichKEGG"))

SelectionWP <- function(RownamesCML, markerList)
  return(MarkersEnrich(markerList, length(RownamesCML), "enrichWP"))

SelectionKEGGEnrich <- function(RownamesCML, markerList)
  return(MarkersEnrich(markerList, length(RownamesCML), "enrichKEGG"))

EnrichmentDescriptions <- function(enrichmentList)
  return(lapply(enrichmentList, function(x){
    #Renaming just to use the p-values in the same way logs are used for marker overlap significance calculations without adding an extra column parameter
    res <- data.frame(avg_log2FC = x$p.adjust)
    rownames(res) <- x$Description
    return(res)
  }))

ClusterPairsOverlap <- function(seuratObj, markers){
  pairs <- c()
  pvalues <- c()
  #shared <- c()
  nClusters <- length(levels(seuratObj))
  nGenes <- length(rownames(seuratObj))
  for (i in 1:(nClusters - 1)){
    message(str_c("Comparing cluster markers overlap for cluster: ", i - 1, "."))
    for (j in (i+1):nClusters){
      pairs <- c(pairs, str_c("Cluster ", i - 1," and Cluster ", j - 1))
      pvalues <- c(pvalues, TwoLogs(markers[[i]], markers[[j]], nGenes))
      #shared <- c(shared, paste(intersect(rownames(markers[[i]]), rownames(markers[[j]])), collapse = ", "))
    }
  }
  df <- data.frame(Grouping = pairs, pvalue = pvalues)
  df <- BYCorrectDF(df)
  return(df)
}

GetNPathways <- function(seuratObj)
  return(length(rownames(GenesER(rownames(seuratObj))@result)))

PathwaysOverlapCluster <- function(seuratObj, descriptions, nTerms){
  pairs <- c()
  pvalues <- c()
  nClusters <- length(levels(seuratObj))
  for (i in 1:(nClusters - 1)){
    message(str_c("Comparing GO pathways overlap for cluster: ", i - 1, "."))
    for (j in (i+1):nClusters){
      pairs <- c(pairs, str_c("Cluster ", i - 1," and Cluster ", j - 1))
      pvalues <- c(pvalues, TwoLogs(descriptions[[i]], descriptions[[j]], nTerms, isLog = F))
    }
  }
  df <- data.frame(Grouping = pairs, pvalue = pvalues)
  df <- BYCorrectDF(df)
  return(df)
}

