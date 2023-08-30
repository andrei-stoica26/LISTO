#' @export

SaveEarlies <- function(gamTradeseq, nKnots, datasetName){
  earlies <- lapply(2:nKnots, function(x) {
    message(str_c("Computing earlie: knot ", x - 1, " to ", x, "."))
    return(earlyDETest(gamTradeseq, knots = c(x - 1, x)))
  }
  )
  saveRDS(earlies, str_c("Earlies", datasetName, ".rds"))
}

SaveStartends <- function(gamTradeseq, datasetName){
  knots <- gamTradeseq@metadata$tradeSeq$knots
  nKnots <- length(knots)
  startends <- lapply(2:nKnots, function(x) {
    message(str_c("Computing startend: knot ", x - 1, " to ", x, "."))
    return(startVsEndTest(gamTradeseq, pseudotimeValues = c(knots[x - 1], knots[x]), lineages = T))
  }
  )
  saveRDS(startends, str_c("Startends", datasetName, ".rds"))
}

SaveCustomStartends <- function(gamTradeseq, bands, datasetName){
  nBands <- length(bands)
  startends <- lapply(2:nBands, function(x) {
    message(str_c("Computing startend: band ", x - 1, " to ", x, "."))
    return(startVsEndTest(gamTradeseq, pseudotimeValues = c(bands[x - 1], bands[x])))
  }
  )
  saveRDS(startends, str_c("CustomStartends", datasetName, ".rds"))
}

SmoothersGrob <- function(gamTradeseq, counts, genes, title, yLim = 5.1){
  dev.new(width = 12, height = 12, noRStudioGD = TRUE)
  return(GrobPlot(lapply(genes, function(x)
    plotSmoothers(gamTradeseq, counts, x, lwd = 0.8) + ggtitle(x) + theme(
      plot.title = element_text(hjust = 0.5, size = 7.5),
      text = element_text(size = 7)) +
      scale_x_continuous(breaks = seq(0, 18, 2)) +
      ylim(0, yLim)
  ),
  nRow = 3, topMargin = 5) + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)))
}

CheckEarlies <- function(earlies, markerList, earlyIndex, markerIndex, earlyPVal = 0.05, markerPVal = 0.05){
  return(earlies[[earlyIndex]][intersect(rownames(subset(markerList[[markerIndex]], p_val_adj < markerPVal)),
                                         rownames(subset(earlies[[earlyIndex]], pvalue < earlyPVal))),])
}

FindEarliesGenes <- function(earlies, markerList, index, pvalue)
  return(lapply(earlies, function(x)intersect(rownames(markerList[[index]]),
                                              rownames(subset(x, pvalue < pvalue)))))

FindUniqueEarliesGenes <- function(earlies, namesList, index, pval = 0.05, fc = 0, wald = 0){
  ueg <- lapply(earlies, function(x) intersect(rownames(subset(x, pvalue < pval &
                                                                 fcMedian > fc & waldStat > wald)),
                                               namesList[[index]]))
  for (i in 2:length(ueg))ueg[[i]] <- setdiff(ueg[[i]], Reduce(union, ueg[1:(i-1)]))
  return(ueg)
}

Knots <- function(gamTradeseq)
  return(data.frame(knotPosition = gamTradeseq@metadata$tradeSeq$knots))

RedYellowPlot <- function(curves, counts, genes, title){
  dev.new(width = 10, height = 12, noRStudioGD = TRUE)
  return(GrobPlot(lapply(genes, function(x)
    plotGeneCount(curves, counts, x) + ggtitle(x) + theme(plot.title = element_text(hjust = 0.5, size = 7.5))),
    nRow = 3, topMargin = 5) + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)))
}
