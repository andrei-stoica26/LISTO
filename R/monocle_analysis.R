#' @export

BuildCellDataSet <- function(seuratObj){
  cellDataSet <- as.cell_data_set(seuratObj)
  cellDataSet <- cluster_cells(cds = cellDataSet, reduction_method = "UMAP")
  cellDataSet <- learn_graph(cellDataSet, use_partition = TRUE)
  return (cellDataSet)
}

OrderCells <- function(cellDataSet, stemCells)
  return(order_cells(cellDataSet, reduction_method = "UMAP", root_cells = stemCells))

PlotTrajectory <- function(monocleRes)
  plot_cells(
    cds = monocleRes,
    color_cells_by = "pseudotime",
    show_trajectory_graph = T,
    label_branch_points = T,
    label_leaves = T,
    label_roots = T,
    graph_label_size = 3,
    trajectory_graph_segment_size = 1.1,
  ) + ggtitle("Monocle pseudotime plot") + theme(plot.title = element_text(hjust = 0.5))

GetPseudotimeCorrelations <- function(seuratObj, method = "spearman")
  return(GetStandardCorrelations(seuratObj, seuratObj$pseudotime, method))

StemSignature <- function(originMarkers, pct1Cutoff, geneCutoff){
  temp <- subset(originMarkers, pct.1 > pct1Cutoff)
  temp$new <- sapply(1:length(rownames(temp)), function(i) temp[i, ]$pct.1/temp[i, ]$pct.2)
  temp <- temp[order(temp$new, decreasing = T),]
  return(rownames(temp)[1:geneCutoff])
}

StemCellsClusters <- function(seuratObj, stemCells)
  return(table(seuratObj$seurat_clusters[stemCells]))

AddPseudotime <- function(seuratObj, monocleRes){
  seuratObj$monocle <- pseudotime(monocleRes)
  return (seuratObj)
}

MonocleAndSignaturePlot <- function(monocleRes, seuratObj, stemSignature){
  dev.new(width =  10, height = 12, noRStudioGD = TRUE)
  p1 <- PlotTrajectory(monocleRes)
  p2 <- Plot_Density_Joint_Only(seuratObj, stemSignature, "viridis") + theme(text =
                                                                               element_text(size = 11))
  return(GrobPlot(list(p1, p2), 2))
}

RunMonocle <- function(seuratObj, originMarkers, pct1Cutoff = 0.4, geneCutoff = 8){
  cds <- BuildCellDataSet(seuratObj)
  stemSignature <- StemSignature(originMarkers, pct1Cutoff, geneCutoff)
  stemCells <- FindCellsCoexpressingGenes(seuratObj, stemSignature)
  print(StemCellsClusters(seuratObj, stemCells))
  monocleRes <- OrderCells(cds, stemCells)
  seuratObj <- AddPseudotime(seuratObj, monocleRes)
  print(MonocleAndSignaturePlot(monocleRes, seuratObj, stemSignature))
  return(seuratObj)
}
