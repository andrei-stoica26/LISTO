#' @export

#These functions operate on the cells expressing genes of interest

GetGenePercentages <- function(seuratObj, genes, toSave = T){
  seurats <- sapply(conditions, function(x) GetAssayData(subset(seuratObj, subset = orig.ident == x),
                                                         slot = "data")[genes,])
  counts <- sapply(seurats, function(x) length(x@p))
  genePercentages <- as.data.frame(mapply(function(x, z)
    lapply(genes, function(y) length(subset(x[y,], x[y,] > 0)) * 100 / z), seurats,counts, SIMPLIFY = T))
  genePercentages <- apply(genePercentages,2,as.character)
  if (toSave)row.names(genePercentages) <- genes
  return(genePercentages)
}

#Use these only for few genes; interrogating sparse matrices is slow
FindCellsExpressingGene <- function (seuratObj, gene, fraction = 0){
  expression <- as.matrix(GetAssayData(seuratObj, slot = "data")[gene,])
  return (rownames(as.data.frame(subset(expression, expression > fraction * max(expression)))))
}

FindNCellsExpressingGene <- function(seuratObj, gene)
  length(FindCellsExpressingGene(seuratObj, gene))

FindPCellsExpressingGene <- function(seuratObj, gene)
  100 * length(FindCellsExpressingGene(seuratObj, gene))/length(colnames(seuratObj))


FindCellsCoexpressingGenes <- function(seuratObj, genes, fraction = 0)
  Reduce(intersect, lapply(genes, function(x) FindCellsExpressingGene(seuratObj, x, fraction)))


SummarizeCellSet <- function(seuratObj, column, cellList = colnames(seuratObj))
  dplyr::count(seuratObj@meta.data[cellList,], {{column}})

#This is to be used for many genes
NonzeroCounts <- function(seuratObj, genes = rownames(seuratObj)){
  expression <- as.matrix(GetAssayData(seuratObj, slot = "data"))
  return (as.data.frame(sapply(genes, function(x) length(which (expression[x,] > 0)))))
  colnames(result) <- c("Number of cells expressing gene")
  return (result)
}

NonzeroPercs <- function(seuratObj, genes = rownames(seuratObj), nCells = length(colnames(seuratObj))){
  expression <- as.matrix(GetAssayData(seuratObj, slot = "data"))
  result <- as.data.frame(sapply(genes, function(x) 100 * length(which (expression[x,] > 0))/nCells))
  colnames(result) <- c("Percentage of cells expressing gene")
  return (result)
}

FilterByPerc <- function(percentageDF, cutoff = 20)
  rownames(subset(percentageDF, percentageDF$`Percentage of cells expressing gene` < cutoff))

AddGeomGenes <- function(seuratObj, genes, colName = "prod"){
  seuratObj@meta.data[[colName]] <- Reduce("*", lapply(genes, function(x) seuratObj@assays$SCT@data[x, ])) ^ (1/length(genes))
  return(seuratObj)
}


NormalizedExpression <- function(seuratObj, gene){
  expression <- as.matrix(GetAssayData(seuratObj, slot = "data")[gene,])
  return(sort(unique(expression)/max(expression), decreasing = T)[-1])
}


GeneCellOverlap <- function(seuratObj, gene1, gene2){
  nCells <- length(colnames(seuratObj))
  discreteExp <- sort(unique(c(NormalizedExpression(seuratObj, gene1), NormalizedExpression(seuratObj, gene2))))
  return (sapply(discreteExp, function(x){
    cells1 <- FindCellsExpressingGene(seuratObj, gene1, x)
    cells2 <- FindCellsExpressingGene(seuratObj, gene2, x)
    return (HyperPValFast(cells1, cells2, nCells))
  }))
}

BYGeneCellOverlap <- function(seuratObj, gene1, gene2)
  return(median(BY(GeneCellOverlap(seuratObj, gene1, gene2), 0.05)$Adjusted.pvalues))

GeneCellOverlapCluster <- function(seuratObj, cluster, gene1, gene2)
  return(GeneCellOverlap(subset(seuratObj, seurat_clusters == cluster), gene1, gene2))

BYGeneCellOverlapCluster <- function(seuratObj, cluster, gene1, gene2)
  return(median(BY(GeneCellOverlap(subset(seuratObj, seurat_clusters == cluster), gene1, gene2), 0.05)$Adjusted.pvalues))

NebulosaMeanPlot <- function(seuratObj, genes){
  seuratObj$nebexp <- colMeans(GetExpressionMatrix(seuratObj)[genes, ])
  return(plot_density(seuratObj, "nebexp"))
}
