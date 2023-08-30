#' @export

#Function used for a stemness-as-gradient analysis in cancer cells.

#Ranking groupings by stemness metric, increasingly or decreasingly as needed
StemnessWilcoxOrig <- function(seuratObj, colStr, altHyp){
  #colStr = "activity", "slingshot" or "monocle
  pairs <- c()
  pvalues <- c()
  colLevels <- levels(seuratObj$orig.ident)
  for (i in colLevels){
    message(str_c("Comparing median ", colStr, " scores for condition: ", i, "."))
    for (j in setdiff(colLevels, i)){
      pairs <- c(pairs, str_c("[", i,"] vs. [", j, "]"))
      pvalues <- c(pvalues, wilcox.test(subset(seuratObj, orig.ident == i)@meta.data[[colStr]],
                                        subset(seuratObj, orig.ident == j)@meta.data[[colStr]],
                                        alternative = altHyp)$p.value)
    }
  }
  df <- data.frame(Grouping = pairs, pvalue = pvalues)
  df <- BYCorrectDF(df)
  return(df)
}

ActivityWilcoxOrig <- function(seuratObj)
  return(StemnessWilcoxOrig(seuratObj, "activity", "greater"))

PseudotimeWilcoxOrig <- function(seuratObj, colStr)
  return(StemnessWilcoxOrig(seuratObj, colStr, "less"))

SlingshotWilcoxOrig <- function(seuratObj)
  return(PseudotimeWilcoxOrig(seuratObj, "slingshot"))

MonocleWilcoxOrig <- function(seuratObj)
  return(PseudotimeWilcoxOrig(seuratObj, "monocle"))

StemnessWilcoxOrigWithinClusters <- function(seurats, seuratObj, typeNames, colStr, altHyp){
  pairs <- c()
  pvalues <- c()
  clusters <- c()
  colLevels <- levels(seuratObj$orig.ident)
  for (cluster in 1:length(seurats)){
    message(str_c("Comparing median ", colStr, " scores for cluster: ", typeNames[cluster], "."))
    for (i in colLevels){
      message(str_c("Comparing median ", colStr, " scores for condition: ", i, "."))
      for (j in setdiff(colLevels, i)){
        clusters <- c(clusters, typeNames[cluster])
        pairs <- c(pairs, str_c("[", i,"] vs. [", j, "]"))
        pvalues <- c(pvalues, wilcox.test(subset(seurats[[cluster]], orig.ident == i)@meta.data[[colStr]],
                                          subset(seurats[[cluster]], orig.ident == j)@meta.data[[colStr]],
                                          alternative = altHyp)$p.value)
      }
    }
  }
  df <- data.frame(Cluster = clusters, Grouping = pairs, pvalue = pvalues)
  df <- BYCorrectDF(df)
  return(df)
}

ActivityWilcoxOrigWithinClusters <- function(seurats, seuratObj, typeNames)
  return(StemnessWilcoxOrigWithinClusters(seurats, seuratObj, typeNames, "activity", "greater"))

SlingshotWilcoxOrigWithinClusters <- function(seurats, seuratObj, typeNames)
  return(StemnessWilcoxOrigWithinClusters(seurats, seuratObj, typeNames, "slingshot", "less"))

MonocleWilcoxOrigWithinClusters <- function(seurats, seuratObj, typeNames)
  return(StemnessWilcoxOrigWithinClusters(seurats, seuratObj, typeNames, "monocle", "less"))

RemoveNullGenes <- function(seuratObj){
  nullGenes <- FindRareGenes(seuratObj, "SCT", 1)
  seuratObj <- RemoveRareGenes(seuratObj, nullGenes)
  return(seuratObj)
}

#Finding max, median, mean and min of stemness metric for each gene, over each cell expressing the gene.

StemnessStats <- function(seuratObj, colStr, ordStat, isDecreasing, doRemoveNull = F){
  x <- Sys.time()
  if (doRemoveNull) seuratObj <- RemoveNullGenes(seuratObj) #Needed when applying it to subsets
  maxp <- c()
  medianp <- c()
  meanp <- c()
  minp <- c()
  ncells <- c()
  for (gene in rownames(seuratObj)){
    cells <- FindCellsExpressingGene(seuratObj, gene)
    values <- seuratObj[[colStr]][cells, ]
    maxp <- c(maxp, max(values))
    medianp <- c(medianp, median(values))
    meanp <- c(meanp, mean(values))
    minp <- c(minp, min(values))
    ncells <- c(ncells, length(cells))
  }
  df <- data.frame(maxp, medianp, meanp, minp, ncells)
  editedColStr <- str_to_title(colStr)
  colnames(df) <- c(str_c(c("Max", "Median", "Mean", "Min"), editedColStr), "NCells")
  rownames(df) <- rownames(seuratObj)
  df <- df[order(df[[str_c(ordStat, editedColStr)]], decreasing = isDecreasing),]
  y <- Sys.time()
  print(y - x)
  return (df)
}

ActivityStats <- function(seuratObj, doRemoveNull = F)
  return (StemnessStats(seuratObj, "activity", "Min", T))

PseudotimeStats <- function(seuratObj, colStr, doRemoveNull = F)
  return (StemnessStats(seuratObj, colStr, "Max", F))

SlingshotStats <- function(seuratObj, doRemoveNull = F)
  return(PseudotimeStats(seuratObj, "slingshot"))

MonocleStats <- function(seuratObj, doRemoveNull = F)
  return(PseudotimeStats(seuratObj, "monocle"))

SplitWithDroppingZeros <- function(seuratObj)
  return(lapply(levels(seuratObj$orig.ident), function(x) {
    seuratSubset <- subset(seuratObj, orig.ident == x)
    rareGenes <- FindRareGenes(seuratSubset, assay = "SCT", 1)
    seuratSubset <- RemoveRareGenes(seuratSubset, rareGenes)
  }))

StemnessConditionStats <- function(seuratObj, FunStemnessStats){
  seurats <- SplitWithDroppingZeros(seuratObj)
  return(lapply(seurats, FunStemnessStats))
}

ActivityConditionStats <- function(seuratObj)
  return(StemnessConditionStats(seuratObj, ActivityStats))

SlingshotConditionStats <- function(seuratObj)
  return(StemnessConditionStats(seuratObj, SlingshotStats))

MonocleConditionStats <- function(seuratObj)
  return(StemnessConditionStats(seuratObj, MonocleStats))

#Plotting number of cells expressing each gene, grouped by the min or max of the stemness-linked metric
StemnessExtremaPlot <- function(stemnessStats, statCol, statStr)
  return(ggplot(stemnessStats, aes(x = {{statCol}}, y = NCells)) +
           geom_point(color =  adjustcolor("goldenrod1", alpha.f = 0.2)) + theme_classic() +
           theme(legend.position = 'None', text = element_text(size = 7), plot.title = element_text(hjust = 0.5)) +
           ggtitle(str_c("Number of cells expressing genes grouped by ", statStr)) +
           labs(y = "Number of cells expressing gene", x = simpleCap(statStr))
  )

ActivityExtremaPlot <- function(stemnessStats, statCol, lowerStatStr)
  return(StemnessExtremaPlot(stemnessStats, {{statCol}}, str_c(lowerStatStr, " ORIGINS activity")))
ActivityMaxPlot <- function(stemnessStats)
  return(ActivityExtremaPlot(stemnessStats, MaxActivity, "maximum"))
ActivityMinPlot <- function(stemnessStats)
  return(ActivityExtremaPlot(stemnessStats, MinActivity, "minimum"))

SlingshotExtremaPlot <- function(stemnessStats, statCol, lowerStatStr)
  return(StemnessExtremaPlot(stemnessStats, {{statCol}}, str_c(lowerStatStr, " Slingshot pseudotime")))
SlingshotMaxPlot <- function(stemnessStats)
  return(SlingshotExtremaPlot(stemnessStats, MaxSlingshot, "maximum"))
SlingshotMinPlot <- function(stemnessStats)
  return(SlingshotExtremaPlot(stemnessStats, MinSlingshot, "minimum"))

MonocleExtremaPlot <- function(stemnessStats, statCol, lowerStatStr)
  return(StemnessExtremaPlot(stemnessStats, {{statCol}}, str_c(lowerStatStr, " Monocle pseudotime")))
MonocleMaxPlot <- function(stemnessStats)
  return(MonocleExtremaPlot(stemnessStats, MaxMonocle, "maximum"))
MonocleMinPlot <- function(stemnessStats)
  return(MonocleExtremaPlot(stemnessStats, MinMonocle, "minimum"))

StemnessExtremaGrob <- function(activityStats, slingshotStats){
  dev.new(width =  8.5, height = 12, noRStudioGD = TRUE)
  plots <- list(ActivityMinPlot(activityStats),
                SlingshotMaxPlot(slingshotStats))
  return(GrobPlot(plots, 2, bottomMargin = 5))
}

#Finding characteristic gene sets (local peaks)

GreedyMat <- function(stemnessStats, stat, stemnessType, isDecreasing, minCells = 300){
  statStr <- str_c(stat, stemnessType)
  stemnessStats <- stemnessStats[order(stemnessStats[[statStr]], decreasing = isDecreasing),]
  genes <- c()
  ncells <- c()
  stats <- c()
  positions <- c()
  target <- minCells
  for (i in 1:nrow(stemnessStats)){
    if (stemnessStats[i,]$NCells > target){
      target <- stemnessStats[i,]$NCells + 1
      genes <- c(genes, rownames(stemnessStats)[i])
      stats <- c(stats, stemnessStats[i,][[statStr]])
      ncells <- c(ncells, target - 1)
      positions <- c(positions, i)
    }
  }
  df <- data.frame(stats, ncells, positions)
  colnames(df) <- c(statStr, "NCells", "Position")
  rownames(df) <- genes
  return(df)
}


ActivityMat <- function(stemnessStats, stat, minCells = 300)
  return(GreedyMat(stemnessStats, stat, "Activity", T, minCells))

PseudotimeMat <- function(stemnessStats, stat, stemnessType, minCells = 300)
  return(GreedyMat(stemnessStats, stat, stemnessType, F, minCells))

SlingshotMat <- function(stemnessStats, stat, minCells = 300)
  return(PseudotimeMat(stemnessStats, stat, "Slingshot", minCells = 300))

MonocleMat <- function(stemnessStats, stat, minCells = 300)
  return(PseudotimeMat(stemnessStats, stat, "Monocle", minCells = 300))

#Plotting characteristic gene sets
LocalPeaksPlot <- function(stemnessStats, statCol, greedyMat, statStr)
  return(ggplot(stemnessStats, aes(x = {{statCol}}, y = NCells)) +
           geom_point(color =  adjustcolor("goldenrod1", alpha.f = 0.2)) + theme_classic() +
           geom_text_repel(data = greedyMat, aes(label = rownames(greedyMat)), size = 2, max.overlaps = 40) +
           geom_point(data = greedyMat,
                      aes(x = {{statCol}}, y = NCells, color = "red"), size = 1.4) +
           geom_smooth(data = greedyMat, se = FALSE, method = "gam", formula = y ~ s(log(x))) +
           theme(legend.position = 'None', text = element_text(size = 7), plot.title = element_text(hjust = 0.5)) +
           ggtitle(str_c("Number of cells expressing genes grouped by ", statStr)) +
           labs(y = "Number of cells expressing gene", x = simpleCap(statStr))
  )

ActivityLocalPeaksPlot <- function(stemnessStats, statCol, lowerStatStr)
  return(LocalPeaksPlot(stemnessStats, {{statCol}},
                        ActivityMat(stemnessStats, simpleCap(lowerStatStr)),
                        str_c(lowerStatStr, " ORIGINS activity")))

ActivityMedianLPP <- function(stemnessStats)
  return(ActivityLocalPeaksPlot(stemnessStats, MedianActivity, "median"))

ActivityMeanLPP <- function(stemnessStats)
  return(ActivityLocalPeaksPlot(stemnessStats, MeanActivity, "mean"))

SlingshotMeanLPP <- function(stemnessStats)
  return(SlingshotLocalPeaksPlot(stemnessStats, MeanSlingshot, "mean"))

SlingshotLocalPeaksPlot <- function(stemnessStats, statCol, lowerStatStr)
  return(LocalPeaksPlot(stemnessStats, {{statCol}},
                        SlingshotMat(stemnessStats, simpleCap(lowerStatStr)),
                        str_c(lowerStatStr, " Slingshot pseudotime")))


SlingshotMedianLPP <- function(stemnessStats)
  return(SlingshotLocalPeaksPlot(stemnessStats, MedianSlingshot, "median"))

SlingshotMeanLPP <- function(stemnessStats)
  return(SlingshotLocalPeaksPlot(stemnessStats, MeanSlingshot, "mean"))

MonocleLocalPeaksPlot <- function(stemnessStats, statCol, lowerStatStr)
  return(LocalPeaksPlot(stemnessStats, {{statCol}},
                        MonocleMat(stemnessStats, simpleCap(lowerStatStr)),
                        str_c(lowerStatStr, " Monocle 3 pseudotime")))

MonocleMedianLPP <- function(stemnessStats)
  return(MonocleLocalPeaksPlot(stemnessStats, MedianMonocle, "median"))

MonocleMeanLPP <- function(stemnessStats)
  return(MonocleLocalPeaksPlot(stemnessStats, MeanMonocle, "mean"))

LocalPeaksGrob <- function(stemnessStats, MedianLPP, MeanLPP){
  dev.new(width =  8.5, height = 12, noRStudioGD = TRUE)
  return(GrobPlot(list(MedianLPP(stemnessStats), MeanLPP(stemnessStats)), 2, bottomMargin = 5))
}

ActivityLPG <- function(stemnessStats)
  return(LocalPeaksGrob(stemnessStats, ActivityMedianLPP, ActivityMeanLPP))

SlingshotLPG <- function(stemnessStats)
  return(LocalPeaksGrob(stemnessStats, SlingshotMedianLPP, SlingshotMeanLPP))

MonocleLPG <- function(stemnessStats)
  return(LocalPeaksGrob(stemnessStats, MonocleMedianLPP, MonocleMeanLPP))

#Finding local peak genes for activity and slingshot
FindLocalPeakGenes <- function(activityStats, slingshotStats)
  return(lapply(list(ActivityMat(activityStats, "Median"),
                     SlingshotMat(slingshotStats, "Median"),
                     ActivityMat(activityStats, "Mean"),
                     SlingshotMat(slingshotStats, "Mean")), rownames))

#Finding out how many times each gene appears among the characteristic gene sets
LocalPeaksFrequencyTable <- function(lpGenes, freqCutoff = 2){
  lpFreqTable <- data.frame(sort(table(unlist(lpGenes)), decreasing = T))
  colnames(lpFreqTable) <- c("Gene", "Frequency")
  lpFreqTable <- subset(lpFreqTable, Frequency > freqCutoff)
  return(lpFreqTable)
}

#Showing intersetions of characteristic gene sets
StatLPUpset <- function(lpGenes){
  dev.new(width =  10, height = 8, noRStudioGD = TRUE)
  allLPGenes <- Reduce(union, lpGenes)
  upsetInput <- as.data.frame(lapply(lpGenes, function(x) allLPGenes %in%x))
  colnames(upsetInput) <- str_c(c("Activity median", "Slingshot pseudotime median",
                                  "Activity mean", "Slingshot pseudotime mean"), " peak genes")
  rownames(upsetInput) <- allLPGenes
  return(UpsetPlot(upsetInput,
                   str_c("Local peak genes sets"),
                   str_c("Intersections of local peak gene sets"),
                   maxDegree = 4
  ))
}

#Computing significance of overlaps of characteristic gene sets
LPGenesOverlap <- function(seuratObj, lpGenes){
  statNames <- c("Activity medians set", "Slingshot pseudotime medians set",
                 "Activity means set", "Slingshot pseudotime means set")
  pairs <- c()
  pvalues <- c()
  for (i in 1:3)
    for (j in (i+1):4){
      pairs <- c(pairs, str_c(statNames[i], " and ", statNames[j]))
      pvalues <- c(pvalues, TwoGeneSetsPVal(lpGenes[[i]], lpGenes[[j]], nrow(seuratObj)))
    }
  df <- data.frame(Grouping = pairs, pvalue = pvalues)
  return (BYCorrectDF(df))
}

GetPeakSetLabels <- function(assays = c("Activity ", "Slingshot pseudotime "),
                             stats = c("medians set", "means set"))
  return(unlist(lapply(stats, function(x) unlist(lapply(assays, function(y)str_c(y, x))))))


#Cut-off for finding high stemness cells. Allowed to vary between top 0.2% and top 16%.

StemnessFilters <- function(seuratObj, colStr, isDecreasing, lowFilter = 0.2, highFilter = 16){
  nCells <- length(colnames(seuratObj))
  return(unique(sort(seuratObj@meta.data[[colStr]],
                     decreasing = isDecreasing)[as.integer(nCells * lowFilter/100):as.integer(nCells * highFilter/100)]))
}

ActivityFilters <- function(seuratObj, lowFilter = 0.2, highFilter = 16)
  return(StemnessFilters(seuratObj, "activity", isDecreasing = T, lowFilter, highFilter))

SlingshotFilters <- function(seuratObj, lowFilter = 0.2, highFilter = 16)
  return(StemnessFilters(seuratObj, "slingshot", isDecreasing = F, lowFilter, highFilter))

MonocleFilters <- function(seuratObj, lowFilter = 0.2, highFilter = 16)
  return(StemnessFilters(seuratObj, "monocle", isDecreasing = F, lowFilter, highFilter))

ColumnCells <- function(seuratObj, column)
  return(dplyr::count(seuratObj@meta.data, {{column}})$n)

IsGreater <- function(x, y)
  return (x > y)

IsLess <- function(x, y)
  return (x < y)

StemnessPValues <- function(seuratObj, index, FiltersFun, stemnessColStr, groupingCol, columnCells,
                            lowerTail, CompFun, lowFilter = 0.2, highFilter = 16){
  stemnessFilters <- FiltersFun(seuratObj, lowFilter, highFilter)
  return(unlist(lapply(stemnessFilters, function(y){
    stemnessCells <- dplyr::count(subset(seuratObj@meta.data, CompFun(seuratObj@meta.data[, stemnessColStr], y)),
                                  {{groupingCol}}, .drop = FALSE)$n
    stemnessTotal <- sum(stemnessCells)
    return(phyper(stemnessCells[index] - (1 - lowerTail), columnCells[index],
                  sum(columnCells) - columnCells[[index]], stemnessTotal, lower.tail = lowerTail))
  })))
}

ActivityPValues <- function(seuratObj, index, groupingCol, columnCells, lowerTail, lowFilter = 0.2, highFilter = 16)
  return(StemnessPValues(seuratObj, index, ActivityFilters, "activity", {{groupingCol}}, columnCells,
                         lowerTail, IsGreater, lowFilter = 0.2, highFilter = 16))

SlingshotPValues <- function(seuratObj, index, groupingCol, columnCells, lowerTail, lowFilter = 0.2, highFilter = 16)
  return(StemnessPValues(seuratObj, index, SlingshotFilters, "slingshot", {{groupingCol}}, columnCells,
                         lowerTail, IsLess, lowFilter = 0.2, highFilter = 16))

MonoclePValues <- function(seuratObj, index, groupingCol, columnCells, lowerTail, lowFilter = 0.2, highFilter = 16)
  return(StemnessPValues(seuratObj, index, MonocleFilters, "monocle", {{groupingCol}}, columnCells,
                         lowerTail, IsLess, lowFilter = 0.2, highFilter = 16))

TopStemnessGroupings <- function(seuratObj, PValuesFun, groupingCol, groupingNames, lowerTail){
  #lowerTail = FALSE: Assessing overrepresentation
  #lowerTail = TRUE: Assessing underrepresentation
  columnCells <- ColumnCells(seuratObj, {{groupingCol}})
  df <- data.frame(Grouping = groupingNames, pvalue = unlist(lapply(1:length(groupingNames), function(x) median(BY(
    PValuesFun(seuratObj, x, {{groupingCol}}, columnCells, lowerTail), 0.05)$Adjusted.pvalues))))
  return(BYCorrectDF(df))
}

TopActivityGroupings <- function(seuratObj, groupingCol, groupingNames, lowerTail)
  return(TopStemnessGroupings(seuratObj, ActivityPValues, {{groupingCol}}, groupingNames, lowerTail))

TopSlingshotGroupings <- function(seuratObj, groupingCol, groupingNames, lowerTail)
  return(TopStemnessGroupings(seuratObj, SlingshotPValues, {{groupingCol}}, groupingNames, lowerTail))

TopMonocleGroupings <- function(seuratObj, groupingCol, groupingNames, lowerTail)
  return(TopStemnessGroupings(seuratObj, MonoclePValues, {{groupingCol}}, groupingNames, lowerTail))
