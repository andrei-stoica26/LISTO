#' @export

#Comparing markers differentially expressed between two clusters
SymmMarkers <- function(seuratObj, cluster1, cluster2){
  nClusters <- length(levels(seuratObj))
  nTests <- nClusters * (nClusters - 1)
  print(nTests)
  markersAB <- FindMarkers(seuratObj, ident.1 = cluster1, ident.2 = cluster2, only.pos = T, min.pct = 0,
                           logfc.threshold = 0, densify = T)
  markersAB <- BonferroniAndFilter(markersAB, nTests)
  markersBA <- FindMarkers(seuratObj, ident.1 = cluster2, ident.2 = cluster1, only.pos = T, min.pct = 0,
                           logfc.threshold = 0, densify = T)
  markersBA <- BonferroniAndFilter(markersBA, nTests)
  return(list(markersAB, markersBA))
}

CellsExpressingSymmMarkers <- function(seuratObj, symmMarkers, cluster1, cluster2){
  seuratSub <- subset(seuratObj, seurat_clusters %in% c(cluster1, cluster2))
  symmNames <- lapply(symmMarkers, rownames)
  return(list(
    sapply(symmNames[[1]], function(gene) FindCellsExpressingGene(seuratSub, gene)),
    sapply(symmNames[[2]], function(gene) FindCellsExpressingGene(seuratSub, gene))
  ))
}

SymmCellsOverlap <- function(seuratObj, symmCells, cluster1, cluster2){
  genes1 <- c()
  genes2 <- c()
  pvalues <- c()
  names1 <- names(symmCells[[1]])
  names2 <- names(symmCells[[2]])
  nCells <- dim(subset(seuratObj, seurat_clusters %in% c(cluster1, cluster2)))[2]
  for (i in 1:length(names1))
    for (j in 1:length(names2)){
      genes1 <- c(genes1, names1[i])
      genes2 <- c(genes2, names2[j])
      pvalues <- c(pvalues, TwoCellSetsPVal(symmCells[[1]][[i]], symmCells[[2]][[j]], nCells))
    }
  df <- data.frame(Gene1 = genes1, Gene2 = genes2, pvalue = pvalues)
  df <- BYCorrectDF(df)
  return(df)
}

OrderedFreqTable <- function(df){
  tempTable <- data.frame(table(df))
  df <- data.frame(Frequency = tempTable$Freq)
  rownames(df) <- tempTable[, 1]
  df <- df[order(df$Frequency, decreasing = T),,drop = F]
  return(df)
}

SymmWrapper <- function(seuratObj, cluster1, cluster2){
  symmMarkers <- SymmMarkers(seuratObj, cluster1, cluster2)
  symmCells <- CellsExpressingSymmMarkers(seuratObj, symmMarkers, cluster1, cluster2)
  symmCellsOverlap <- SymmCellsOverlap(seuratObj, symmCells, cluster1, cluster2)
  symmDF1 <- OrderedFreqTable(symmCellsOverlap$Gene1)
  symmDF2 <- OrderedFreqTable(symmCellsOverlap$Gene2)
  return(list(symmDF1, symmDF2))
}
