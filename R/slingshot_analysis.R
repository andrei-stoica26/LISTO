#' @export

#Slingshot functions
GetDimred <- function(seuratObj)
  return(seuratObj@reductions$umap@cell.embeddings)

GetSlingshotCurves <- function(seuratObj, dimred){
  lineages <- getLineages(data = dimred, clusterLabels = seuratObj$seurat_clusters)
  curves <- getCurves(lineages)
  return (curves)
}

PlotLineages <- function(seuratObj, curves, dimred){
  plot(dimred, col = hue_pal()(length(levels(seuratObj)))[seuratObj$seurat_clusters], cex = 0.1, pch = 16)
  lines(SlingshotDataSet(curves))
  p <- GrabGrob()
  p <- as.ggplot(p)
  return(p + ggtitle("Differentiation lineages") + theme(plot.title = element_text(hjust = 0.5, vjust = -12)))
}

AddSlingshotResults <- function(seuratObj, curves, colStr, fun){
  nLineages <- ncol(fun(curves))
  for (i in 1:nLineages)
    #replace_na to prevent issues with VlnPlot - if the first cell is NA then VlnPlot will give an error
    seuratObj@meta.data[[str_c(colStr, i)]] <- replace_na(fun(curves)[, i], NaN)
  return (seuratObj)
}

AddLineages <- function(seuratObj, curves)
  return(AddSlingshotResults(seuratObj, curves, "Lineage", slingPseudotime))

AddCurveweights <- function(seuratObj, curves)
  return(AddSlingshotResults(seuratObj, curves, "Curveweight", slingCurveWeights))


SaveLineageGam <- function(seuratLineage, colStr){
  expressionLineage <- GetExpressionMatrix(seuratLineage)
  rareGenes <- FindRareGenes(seuratLineage, "SCT", 10)
  seuratLineage <- RemoveRareGenes(seuratLineage, rareGenes)
  x <- Sys.time()
  gam <- FitGam(seuratLineage, expressionLineage, rownames(seuratLineage), colStr)
  y <- Sys.time()
  print(y - x)
  SaveGam(gam, str_c(datasetName, colStr, "GAM.rds"))
}

#Useful for adding labels on a figure
FindClusterCenters <- function(seuratObj, dimred){
  x <- c()
  y <- c()
  for (i in levels(seuratObj)){
    clusterCells <- colnames(subset(seuratObj, seurat_clusters == i))
    x <- c(x, mean(dimred[clusterCells, 1]))
    y <- c(y, mean(dimred[clusterCells, 2]))
  }
  return(data.frame(x, y))
}

LineagePlot <- function(seuratObj, dimred, curves, title, clCenters){
  lineage <- gsub(" ", "", title)
  dev.new(width =  7, height = 7, noRStudioGD = TRUE)
  p1 <- DimPlot(seuratObj, label = T) + theme(text = element_text(size = 7),
                                              plot.title = element_text(hjust = 0.5)) + ggtitle("Clusters")
  #Filtering out NaN with the > -1 comparison
  expr <- FetchData(seuratObj, vars = lineage)
  p2 <- VlnPlot(seuratObj[, which (expr > - 1)], lineage) + theme(text = element_text(size = 7)) +
    ggtitle (str_c("Violin plot - ", title, " pseudotime"))
  p3 <- FeaturePlot(seuratObj, lineage) + theme(text = element_text(size = 7)) +
    ggtitle (str_c("Feature plot - ", title, " pseudotime"))
  print(length(unique(hue_pal()(length(levels(seuratObj)))[seuratObj$seurat_clusters])))
  plot(dimred, col = hue_pal()(length(levels(seuratObj)))[seuratObj$seurat_clusters], cex = 0.1, pch = 16) +
    text(x = clCenters$x, y = clCenters$y, labels = levels(seuratObj), cex = 0.8)
  t <- SlingshotDataSet(curves)
  t@lineages <- list(t@lineages[[lineage]])
  t@curves <- list(t@curves[[lineage]])
  lines(t)
  #Axis labels disappear in p4, axis.title.x/y is needed to bring them back
  p4 <- as.ggplot(GrabGrob()) +
    theme(text = element_text(size = 7), axis.title.x = element_text(),
          axis.title.y = element_text(angle = 90, vjust = 5),
          plot.title = element_text(face = "bold", hjust = 0.5, vjust = -8)) +
    labs(x = "UMAP_1", y = "UMAP_2") + ggtitle ("Differentiation trajectory")
  return (GrobPlot(list(p1, p2, p3, p4), nRow = 2, c(20, 15, 5, 15)) + ggtitle(str_c(title, " ordering")) +
            theme(plot.title = element_text(hjust = 0.5)))
}

GenePseudotime <- function(seuratObj, curves, gene)
  return(slingPseudotime(curves)[FindCellsExpressingGene(seuratObj, gene),])

FilterPseudotime <- function(pt, lineage, high = 1, low = 0){
  filteredCells <- c()
  for (cell in rownames(pt))
    if (!is.na(pt[cell, lineage]) & pt[cell, lineage] < high &
        pt[cell, lineage] >= low)
      filteredCells <- c(filteredCells, cell)
  return(pt[filteredCells, ])
}

FilterPseudotimeWrapper <- function(seuratObj, curves, gene, lineage, high = 1, low = 0)
  return(FilterPseudotime(GenePseudotime(seuratObj, curves, gene), lineage, high, low))

LowPseudotimeCells <- function(curves, high){
  pt <- slingPseudotime(curves)
  lineages <- colnames(pt)
  cells <- c()
  for (lineage in lineages)
    cells <- c(cells, rownames(FilterPseudotime(pt, lineage, high)))
  return(unique(cells))
}

MaxGenePseudotime <- function(gp){
  gp[is.na(gp)] <- -1
  return(max(gp))
}

MinGenePseudotime <- function(gp){
  gp[is.na(gp)] <- 9999
  return(min(gp))
}

MaxPseudotimeAllGenes <- function(seuratObj, curves)
{
  x <- Sys.time()
  mp <- c()
  ncells <- c()
  for (gene in rownames(seuratObj)){
    currentCells <- FindCellsExpressingGene(seuratObj, gene)
    mp <- c(mp, MaxGenePseudotime(slingPseudotime(curves)[currentCells, ]))
    ncells <- c(ncells, length(currentCells))
  }
  df <- data.frame(MaxPseudotime = mp, NCells = ncells)
  rownames(df) <- rownames(seuratObj)
  df <- df[order(df$MaxPseudotime),]
  y <- Sys.time()
  print(y - x)
  return (df)
}

ProdNA <- function(x, y)
  return (if (NA %in% c(x, y)) 0 else x * y)

WeightedMeanNA <- function(values, weights){
  #In the input, values can contain NA; weights cannot
  return(sum(sapply(1:length(values), function(i) ProdNA(values[[i]], weights[[i]])))/sum(weights))
}

AddCombinedPseudotime <- function(seuratObj, curves){
  seuratObj$slingshot <- sapply(rownames(slingPseudotime(curves)), function(cell)
    WeightedMeanNA(slingPseudotime(curves)[cell, ], slingCurveWeights(curves)[cell, ]))
  return (seuratObj)
}

