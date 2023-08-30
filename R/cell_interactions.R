#' @export

GetCellChat <- function(seuratObj){
  data.input <- GetAssayData(seuratObj, assay = "SCT", slot = "data")
  meta <- data.frame(group = seuratObj$cell.type, row.names = colnames(seuratObj))
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  cellchat@DB <- CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  cellchat@var.features[["features"]] <- VariableFeatures(seuratObj)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = F, type = "truncatedMean", trim = 0.1)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}

NetVisualAgg <- function(cellchat, pathways){
  dev.new(width =  12, height = 12, noRStudioGD = TRUE)
  return(GrobPlot(lapply(pathways, function(x){
    netVisual_aggregate(cellchat, x, vertex.label.cex = 0.6, pt.title = 2, edge.width.max = 2,
                        title.space = 10000, vertex.size.max = 2)
    p <- GrabGrob()
    return(as.ggplot(p))
  }), nRow = 4, topMargin = 5) +
    ggtitle("Signalling pathways found between clusters") + theme(plot.title = element_text(hjust = 0.5)))
}

