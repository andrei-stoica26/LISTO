#' @export

AddNumToList <- function (num_list, val) lapply(num_list, function(x) x + val)

SetIdents <- function (ident, type, conds = NULL){
  if (type == "Condition"){
    ident <- lapply(ident, function(ident) match(ident, conds))
    return (AddNumToList (ident, 1))
  }
  return (AddNumToList (ident, 2))
}

SetReactomeColumnNames <- function(reactome, type){
  if (type == "Cluster")
    colnames(reactome@results$Seurat$pathways) <- gsub("X", "Cluster ", colnames(reactome@results$Seurat$pathways))
  return (reactome)
}

CreateReactome <- function(seuratObj, type){
  seuratObj <- RenameAssays(seuratObj, SCT = "RNA")
  Idents(seuratObj) <- type
  reactome <- analyse_sc_clusters(seuratObj, verbose = TRUE)
  reactome <- SetReactomeColumnNames(reactome, type)
  return (reactome)
}

setClass ("ExtremaFunctions",
          slots=list(colSel = "function", best1 = "function", best2 = "function", sign = "numeric"))

min.col <- function(m) max.col(-m)

AddIdentColumns <- function (dpe, ident.str, ident.vector, ext.fxn){
  if (length(ident.vector) > 1)
    dpe[, ident.str] <- apply(dpe[, -1][, (ident.vector - 1)], 1, ext.fxn)
  else
    dpe[, ident.str] <- dpe[, ident.vector]
  return(dpe)
}

DEPathways <- function(pe, id1.vector, id2.vector, direction){
  if (direction == "up")
    EF <- new("ExtremaFunctions", colSel = max.col, best1 = min, best2 = max, sign = 1)
  else
    EF <- new("ExtremaFunctions", colSel = min.col, best1 = max, best2 = min, sign = -1)

  # Save only the pathways where the expression level reaches the desired extreme point in the selection
  dpe <- subset(pe, EF@colSel(pe[, -1]) %in% (id1.vector - 1))


  dpe <- AddIdentColumns(dpe, "ident.1", id1.vector, EF@best1)
  dpe <- AddIdentColumns(dpe, "ident.2", id2.vector, EF@best2)
  dpe[, "diff"] <- EF@sign * (dpe$ident.1 - dpe$ident.2)
  dpe <- arrange(dpe[, c("Name", "ident.1", "ident.2", "diff")], -diff)
  return(dpe)
}

PathwayImageName <- function(type, id1, id2, direction){
  return (str_c(direction, "regulated in ", ifelse (type == "Cluster", str_c("cluster ", toString(id1 - 2), " vs. all other clusters"),
                                                    str_c(paste(unlist(conditions[id1 - 1]), collapse=' & '), " vs. ",
                                                          paste(unlist(conditions[id2 - 1]), collapse=' & ') ))))
}

PathwayPlot <- function(reactome, dpe,id1, id2, type, direction, n_colors){
  caption <- str_c("Pathways ", PathwayImageName(type, id1, id2, direction), ".")
  numbers <- AddFileNumber(1)
  plots <- lapply(rownames(dpe)[1:4], function(x) plot_gsva_pathway(reactome, pathway_id = x) +
                    geom_bar(stat="identity",
                             color = brewer.pal(n_colors, "Set3"),
                             fill = brewer.pal(n_colors, "Set3")) +
                    ylim(-1, 1) +
                    xlab(type) +
                    theme(plot.title = element_text(size = 13 - 0.3 * n_colors))
  )
  PDFSlide(str_c(caption, "pdf"), GrobPlot(plots, 2, caption), numbers)
}

PathwayHeatmapPlot <- function(reactome, type){
  plot_gsva_heatmap(reactome, margins = c(1.5, 10), truncate_names = F, cexCol = 0.7, cexRow = 0.7, srtCol = 0,
                    adjCol = c(NA, 0))
  p <- GrabGrob()
  sharedString <- str_c(tolower(type), "s.")
  PDFSlide(str_c("Pathway heatmap plot - ", sharedString, "pdf"),
           SimplePlot(as.ggplot(p), str_c("Heatmap of the top 20 enriched pathways - ", sharedString),
                      "Figure"), AddFileNumber(1))
}


PathwayAnalysis <- function(reactome, id1, id2, type){
  pe <- pathways(reactome)
  invisible(mapply(function(x, y){
    lapply(c("up", "down"), function(z){
      dpe <- DEPathways(pe, x, y, z)
      PathwayPlot(reactome, dpe, x, y, type, z, length(pe) - 1)
    })
  }, id1, id2))
  PathwayHeatmapPlot(reactome, type)
}

