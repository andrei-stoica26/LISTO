#' @export

#Check the distribution of cells from each treatment condition in clusters
CheckCondDist <- function(seuratObj, isUnderrep){
  seuratObj$orig.ident <- droplevels(seuratObj$orig.ident)
  conditionCells <- dplyr::count(seuratObj@meta.data, orig.ident)$n
  clusterCells <- dplyr::count(seuratObj@meta.data, seurat_clusters)$n
  nCells <- length(colnames(seuratObj))
  condCells <- lapply(levels(seuratObj$orig.ident), function(x)
    dplyr::count(subset(seuratObj@meta.data, orig.ident == x), seurat_clusters)$n)
  OCC <- Reduce(rbind, lapply(1:length(clusterCells), function(x){
    message(str_c("Calculating representation of condition in cluster ", x - 1, "..."))
    return(data.frame(Cluster = rep(x - 1, length(conditionCells)),
                      Condition = levels(seuratObj$orig.ident),
                      PercInCluster = unlist(lapply(1:length(conditionCells),
                                                    function(y) condCells[[y]][x] / clusterCells[x] * 100)),
                      PercOverall = conditionCells/nCells * 100,
                      pvalue = unlist(lapply(1:length(conditionCells),
                                             function(y) phyper(condCells[[y]][x] - (1 - isUnderrep),
                                                                conditionCells[y],
                                                                nCells - conditionCells[y],
                                                                clusterCells[x],
                                                                lower.tail = isUnderrep)))))
  }))
  OCC <- OCC[order(OCC$pvalue), ]
  OCC$pvalue <- BY(OCC$pvalue, 0.05)$Adjusted.pvalues
  OCC <- subset(OCC, OCC$pvalue < 0.05)
  return(OCC)
}

PrintOCC <- function(OCC, types){
  OCC$Cluster <- types[OCC$Cluster + 1]
  MessageVector(OCC$Condition)
  MessageVector(OCC$Cluster)
  PrintVector(OCC$pvalue)
  PrintVector(OCC$PercInCluster)
  PrintVector(OCC$PercOverall)
}

OrderedSplit <- function(seuratObj)
  return(lapply(levels(seuratObj), function(x) return(subset(seuratObj, seurat_clusters == x))))

SingleGeneSelectionsInClusters <- function(singleGenes, CCMarkers, RownamesCML, index,  typeNames){
  print(singleGenes[index])
  clusters <- c()
  pvalues <- c()
  selections <- c()
  for (i in 1:length(CCMarkers))
    for (j in 1:length(RownamesCML)){
      if (singleGenes[[index]] %in% rownames(CCMarkers[[i]][[j]])){
        clusters <- c(clusters, typeNames[i])
        selections <- c(selections, RownamesCML[j])
        pvalues <- c(pvalues, CCMarkers[[i]][[j]][singleGenes[[index]],,drop = F]$p_val_adj)
      }
    }
  df <- data.frame(Cluster = clusters, Grouping = selections, pvalue = pvalues)
  if (length(rownames(df)) > 0)
    df <- df[order(df$pvalue),]
  return(df)
  #It is already corrected with Bonferroni - no need for BY
}

LitSelectionsInClusters <- function(LitMarkers, CCMarkers, RownamesCML, index,  typeNames, N){
  clusters <- c()
  pvalues <- c()
  selections <- c()
  for (i in 1:length(CCMarkers))
    for (j in 1:length(RownamesCML)){
      clusters <- c(clusters, typeNames[i])
      selections <- c(selections, RownamesCML[j])
      pvalues <- c(pvalues, OneLitOneLog(LitMarkers[[index]], CCMarkers[[i]][[j]], N))
    }
  df <- data.frame(Cluster = clusters, Grouping = selections, pvalue = pvalues)
  return(BYCorrectDF(df))
}

LogSelectionsInClusters <- function(logMarkers, CCMarkers, RownamesCML, typeNames, N){
  clusters <- c()
  pvalues <- c()
  selections <- c()
  for (i in 1:length(CCMarkers))
    for (j in 1:length(RownamesCML)){
      clusters <- c(clusters, typeNames[i])
      selections <- c(selections, RownamesCML[j])
      pvalues <- c(pvalues, TwoLogs(logMarkers, CCMarkers[[i]][[j]], N))
    }
  df <- data.frame(Cluster = clusters, Grouping = selections, pvalue = pvalues)
  return(BYCorrectDF(df))
}
