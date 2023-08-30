#' @export

#Functions for getting corrlations between the expression of genes

GetStandardCorrelations <- function(seuratObj, standard, method = "spearman"){
  geneIndices <- 1:length(rownames(seuratObj))
  expMatrix <- as.matrix(seuratObj@assays$SCT@data)[geneIndices, ]
  geneNames <- rownames(seuratObj)[geneIndices]
  result <- as.data.frame(sapply(1:length(geneNames), function (x) round(cor(expMatrix[x,], standard, method = method), 6)), row.names = geneNames)
  colnames(result) <- "Correlation"
  return(result[order(-result$Correlation), , drop = F])
}

GetExpressionMatrix <- function(seuratObj)
  as.matrix(seuratObj@assays$SCT@data)

PairCorrelation <- function(seuratObj, gene1, gene2, method = "spearman"){
  #Slow: use GetExpressionCorrelations for many genes
  expressionMatrix <- GetExpressionMatrix(seuratObj)
  return (stats::cor(expressionMatrix[gene1, ], expressionMatrix[gene2, ]))
}

GetIndex <- function(gene, geneNames)
  which (geneNames == gene)

GetAllExpression <- function(seuratObj, geneNames){
  geneMatrix <- GetExpressionMatrix(seuratObj)
  geneIndices <- sapply(geneNames, function(x) GetIndex(x, rownames(geneMatrix)))
  print(geneIndices)
  return(t(geneMatrix[geneIndices, ]))
}

GetExpressionCorrelations <- function(seuratObj, selectedGene, method = "pearson")
  GetStandardCorrelations(seuratObj, seuratObj@assays$SCT@data[selectedGene, ], method)

GetGeneCorrelationMatrix <- function(seuratObj, genes, method = "spearman"){
  expressionMatrix <- GetExpressionMatrix(seuratObj)[genes, ]
  return (cor(t(expressionMatrix)))
}
