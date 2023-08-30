#' @export

SignatureRepresentation <- function(seuratObj, signature){
  seuratObj$orig.ident <- droplevels(seuratObj$orig.ident)
  sigCells <- FindCellsCoexpressingGenes(seuratObj, signature)
  sigSeurat <- subset(seuratObj, cells = sigCells)
  conditionSeurat <- dplyr::count(seuratObj@meta.data, orig.ident)$n
  conditionSig <- dplyr::count(sigSeurat@meta.data, orig.ident, .drop = FALSE)$n
  print(dplyr::count(sigSeurat@meta.data, orig.ident, .drop = FALSE))
  df <- data.frame(Grouping = levels(seuratObj$orig.ident), pvalue = sapply(1:length(levels(seuratObj$orig.ident)),
                                                                            function(i) phyper(conditionSig[i], conditionSeurat[i], sum(conditionSeurat) - conditionSeurat[i],
                                                                                               sum(conditionSig))))
  df <- BYCorrectDF(df)
  return(df)
}

ClusterSignatureRepresentation <- function(seuratObj, signature){
  seuratObj$seurat_clusters <- droplevels(seuratObj$seurat_clusters)
  sigCells <- FindCellsCoexpressingGenes(seuratObj, signature)
  sigSeurat <- subset(seuratObj, cells = sigCells)
  groupingSeurat <- dplyr::count(seuratObj@meta.data, seurat_clusters)$n
  groupingSig <- dplyr::count(sigSeurat@meta.data, seurat_clusters, .drop = FALSE)$n
  print(dplyr::count(sigSeurat@meta.data, seurat_clusters, .drop = FALSE))
  df <- data.frame(Grouping = levels(seuratObj$seurat_clusters),
                   pvalue = sapply(1:length(levels(seuratObj$seurat_clusters)),
                                   function(i) phyper(groupingSig[i], groupingSeurat[i], sum(groupingSeurat) - groupingSeurat[i], sum(groupingSig))))
  df <- BYCorrectDF(df)
  return(df)
}

GroupingGenes <- function(singleGenes, markers){
  pairs <- c()
  pvalues <- c()
  for (x in singleGenes)
    if (x %in% rownames(markers)){
      pairs <- c(pairs, x)
      pvalues <- c(pvalues, markers[x, ]$p_val_adj)
    }
  #Markers are called Grouping just to be consistent with PrintOverlap2
  df <- data.frame(Grouping = pairs, pvalue = pvalues)
  if (length(df$pvalue) > 0)df <- df[order(df$pvalue),]
  return(df)
  #It is already corrected with Bonferroni - no need for BY!
}

AllGroupingsGenes <- function(singleGenes, markerList)
  return(lapply(markerList, function(x) GroupingGenes(singleGenes, x)))

PrintGenesInGroupings <- function(singleGenes, markerList)
  return(invisible(lapply(AllGroupingsGenes(singleGenes, markerList), PrintOverlap)))

SingleGeneSelections <- function(singleGenes, selectionMarkers, RownamesCML, index){
  print(singleGenes[index])
  pairs <- c()
  pvalues <- c()
  for (i in 1:length(RownamesCML))
    if (singleGenes[[index]] %in% rownames(selectionMarkers[[i]])){
      pairs <- c(pairs, RownamesCML[i])
      pvalues <- c(pvalues, selectionMarkers[[i]][singleGenes[index], ]$p_val_adj)
    }
  df <- data.frame(Grouping = pairs, pvalue = pvalues)
  if (dim(df)[1] > 0)df <- df[order(df$pvalue),]
  return(df)
  #It is already corrected with Bonferroni - no need for BY!
}

PrintAllSGSelections <- function(singleGenes, selectionMarkers, RownamesCML)
  invisible(lapply(1:length(singleGenes), function(x)
    PrintOverlap(SingleGeneSelections(singleGenes, selectionMarkers, RownamesCML, x))))
