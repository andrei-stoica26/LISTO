#' @export

#Tools for working with ORIGINS activity

ActivityFP <- function(seuratObj)
  return(FeaturePlot(seuratObj, "activity") + scale_colour_gradientn(colours = wes_palette("Zissou1")))

ActivityHist <- function(seuratObj,
                         histTitle = "The distribution of activity scores in each experimental condition",
                         nBins = 100)
  return (MakeHistogram(seuratObj, "activity", "orig.ident", "Activity", nBins = 100) + ggtitle(histTitle)+
            theme(plot.title = element_text(hjust = 0.5)))

ActivityHistCl <- function(seuratObj, histTitle = "The distribution of activity scores in each cluster", nBins = 100)
  return (MakeHistogram(seuratObj, "activity", "seurat_clusters", "Activity", nBins = 100) + ggtitle(histTitle)+
            theme(plot.title = element_text(hjust = 18), axis.title.x = element_text(vjust = 1.8)))

GetActivityCorrelations <- function(seuratObj, method = "spearman")
  return(GetStandardCorrelations(seuratObj, seuratObj$activity, method))

TopActivityGenes <- function(activityCorrs, cutoff = 0.2)
  return(rownames(subset(activityCorrs, Correlation > cutoff)))

ActsPlotDF <- function(seuratObj, acts, colStr, df)
  return (data.frame(Gene = acts, Value = df[acts, colStr],
                     Cluster = ClusterAvgExp(seuratObj, acts)$rowmax))

#Plotting genes with a high correlation with activity scores, arranged by the cluster of maximum avg. activity scores
ActsPlot <- function(seuratObj, actsPlotDF, geneCutoff = 50)
  return (ggplot(data = actsPlotDF) + aes(y = Value, x = Cluster) + geom_point(color = "red") +
            geom_text_repel(aes(label = Gene), size = 3, max.overlaps = Inf) + theme_minimal() +
            ggtitle(str_c("Top ", geneCutoff,
                          " genes ranked by the Spearman correlation of their expression \nand ORIGINS activity, grouped by the cluster of highest expression")) +
            theme(plot.title = element_text(hjust = 0.5, size = 12)) +
            coord_cartesian(xlim = c(0, length(levels(seuratObj)))) +
            scale_x_discrete(limits=0:length(levels(seuratObj)) - 1) +
            ylab("Mean weighted activity")
  )

WeightedActsPlot <- function(actsPlotDF, seuratObj, cutoff = 0.5)
  return (ActsPlot(actsPlotDF, seuratObj, cutoff))


#Comparing activity scores between clusters
ActivityWilcoxCluster <- function(seuratObj){
  pairs <- c()
  pvalues <- c()
  colLevels <- levels(seuratObj$seurat_clusters)
  for (i in colLevels){
    message(str_c("Comparing median activity scores for cluster: ", i, "."))
    for (j in setdiff(colLevels, i)){
      pairs <- c(pairs, str_c("Cluster ", i," vs. Cluster ", j))
      pvalues <- c(pvalues, wilcox.test(subset(seuratObj, seurat_clusters == i)$activity, subset(seuratObj, seurat_clusters == j)$activity,
                                        alternative = "greater")$p.value)
    }
  }
  df <- data.frame(Grouping = pairs, pvalue = pvalues)
  df <- BYCorrectDF(df)
  return(df)
}

ScoreActivityWilcox <- function(seuratObj, actWilDF){
  df <- data.frame(Cluster = str_c("Cluster ", levels(seuratObj)),
                   Wins = sapply(levels(seuratObj),
                                 function(x)sum(sapply(actWilDF$Grouping,
                                                       function(y) strsplit(y, " vs. ")[[1]][[1]] ==
                                                         str_c("Cluster ", x)))),
                   Losses = sapply(levels(seuratObj),
                                   function(x)sum(sapply(actWilDF$Grouping,
                                                         function(y) strsplit(y, " vs. ")[[1]][[2]] ==
                                                           str_c("Cluster ", x))))
  )
  df$Score <- df$Wins - df$Losses
  return(df[order(df$Score, decreasing = T),])
}

#Find genes associated with the highest activity score relative to their expression
WAActivity <- function(seuratObj){
  expression <- as.matrix(GetAssayData(seuratObj, assay = "SCT", slot = "data"))
  df <- data.frame(WAActivity = unlist(lapply(rownames(seuratObj), function(x){
    geneExp <- expression[x,]
    return (sum(geneExp * seuratObj$activity) / sum(geneExp))
  }
  )), nCells = unlist(lapply(rownames(seuratObj), function(x) length(FindCellsExpressingGene(seuratObj, x)))))
  rownames(df) <- rownames(seuratObj)
  df <- df[order(df$WAActivity, decreasing = T), ,drop = F]
  return(df)
}

ActivityCohen <- function(seuratObj, cond1, cond2){
  sdev <- sd(seuratObj$activity)
  diff <- mean(subset(seuratObj, orig.ident == cond1)$activity) - mean(subset(seuratObj, orig.ident == cond2)$activity)
  return(diff / sdev)
}


