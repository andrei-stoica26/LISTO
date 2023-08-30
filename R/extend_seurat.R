#' @export

#Functions in this file add columns to Seurat metadata

#Grouping a continuous measurement into categories
SetBandLabels <- function(breakPoints){
  brLength <- length(breakPoints)
  midSections <- unlist(lapply(seq(1:(brLength-1)),
                               function(x) str_c("Between ", breakPoints[x], " and ", breakPoints[x+1])))
  return (c(str_c("Below ", breakPoints[1]), midSections, str_c("Above ", breakPoints[brLength])))
}

#Adding said categories to the metadata
AddBand <- function(seuratObj, bandName, columnName, breakPoints){
  ticks <- SetBandLabels(breakPoints)
  seuratObj@meta.data[[bandName]] <-
    cut(seuratObj@meta.data[[columnName]], breaks = c(-Inf, breakPoints, Inf),
        labels = ticks)
  seuratObj@meta.data[[bandName]] <- factor(seuratObj@meta.data[[bandName]],levels = rev(ticks))
  return(seuratObj)
}

#Adding cell types
AnnotateSeurat <- function(seuratObj, colName, types){
  seuratObj@meta.data[[colName]] <- sapply(seuratObj$seurat_clusters, function(x) types[x])
  return (seuratObj)
}

GeneExpression <- function(expression, gene)
  return (expression[gene, ])

GenePresence <- function(expression, gene)
  return (ifelse(expression[gene, ] > 0, 1, 0))

AppendGeneExpression <- function(seuratObj, expression, gene, colName = "geneExp"){
  seuratObj@meta.data[[colName]] <- GeneExpression(expression, gene)
  return(seuratObj)
}

AppendGenePresence <- function(seuratObj, expression, gene, colName = "genePres"){
  seuratObj@meta.data[[colName]] <- GenePresence(expression, gene)
  return(seuratObj)
}
