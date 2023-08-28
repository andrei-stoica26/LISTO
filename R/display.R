#Creates a grob plot with a defined number of rows, where each subplot is labelled by a letter


GrobPlot <- function(plotObjects, nRow, topMargin = 0, rightMargin = 0, bottomMargin = 0, leftMargin = 0)
  return(as.ggplot(arrangeGrob(grobs = mapply(function(x, y) x + labs(tag = y, fontface = 2)
                                              + theme(plot.margin = margin(topMargin, rightMargin,
                                                                           bottomMargin, leftMargin, "pt"),
                                                      plot.tag = element_text(face = "bold", size = 10))
                                              + coord_cartesian(clip = "off"),
                                              plotObjects, LETTERS[1:length(plotObjects)], SIMPLIFY = F),
                               nrow = nRow)))


#Creates a histogram in which sections can be delineated using one or two lines

MakeHistogram <- function(seuratObj, xColStr, fillColStr, xLabel, palette = wes_palette("Darjeeling1"),
                          nBins = 100, vline1 = NULL, vline2 = NULL){
  ggplot(seuratObj@meta.data, aes(x = !!as.name(xColStr), color = !!as.name(fillColStr),
                                  fill = !!as.name(fillColStr))) +
    geom_histogram(bins = nBins, position="identity", alpha = 0.5) +
    labs(x = xLabel, y = "Cell count") +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.key.size = unit(0.5, 'cm'),
          legend.text=element_text(size = 10)) +
    {if (length(unique(factor(seuratObj@meta.data[[fillColStr]]))) < 6)
      scale_fill_manual(values = palette)} +
    {if (length(unique(factor(seuratObj@meta.data[[fillColStr]]))) < 6)
      scale_colour_manual(values = palette)} +
    {if(!is.null(vline1)) geom_vline(aes(xintercept = vline1), color="blue", linetype="dashed", size=1)} +
    {if(!is.null(vline2)) geom_vline(aes(xintercept = vline2), color="blue", linetype="dashed", size=1)}
}

#Creates the plot of eight histograms used to illustrate quality control choices

EightHistograms <- function(seurats, index){
  dev.new(width =  10, height = 12, noRStudioGD = TRUE)
  titles <- c(str_c(c("Frequency of doublets", strings3[-1]), " - ", conditions[index]))
  vlines <- allBounds[[index]]
  plots <- lapply(1:8, function(x)
    MakeHistogram(seurats[[index]], strings1[x], strings2[x], strings3[x],
                  vline1 = vlines[[x]][1], vline2 = vlines[[x]][2]) +
      ggtitle(titles[x]) +
      theme(text = element_text(size = 7.5), legend.position = if(x==1) "right" else "none",
            legend.key.size = unit(0.2, 'cm'),
            legend.text = element_text(size = 8),
            plot.title = element_text(size = 8, hjust = 0.5),
            axis.title.x = element_text(vjust = 1.8)))
  return(GrobPlot(plots, 4, leftMargin = 15, topMargin = 15))
}

#creates violin and feature plots of a list of genes

VlnFeatTitle <- function(genes, inter = "")
  if (inter == "")
    return(str_c("The expression of ", inter, paste(genes[1:length(genes) - 1], collapse = ", "), " and ",
                 genes[length(genes)])) else return (inter)

VlnFeat <- function(seuratObj, genes, inter = ""){
  dev.new(width =  10, height = 12, noRStudioGD = TRUE)
  violins <- lapply(genes, function(x) VlnPlot(seuratObj, x) + ggtitle(str_c("Violin plot - ", x)) +
                      theme(text = element_text(size = 7), plot.title = element_text(hjust = 0.5)) + NoLegend())

  features <- lapply(genes, function(x) FeaturePlot(seuratObj, x) + ggtitle(str_c("Feature plot - ", x)) +
                       theme(text = element_text(size = 7), plot.title = element_text(hjust = 0.5)))
  rightOrder <- unlist(lapply(1:length(genes), function(x) c(x, x + length(genes))))
  return (GrobPlot(c(features, violins)[rightOrder], length(genes), topMargin = 10) +
            ggtitle(VlnFeatTitle(genes, inter)) + theme(plot.title = element_text(hjust = 0.5, size = 10.5)))
}

#creates violin and feature plots for TDPM genes and titles it accordingly

TDPMVlnFeat <- function(seuratObj, genes)
  return (VlnFeat(seuratObj, genes, str_c("The expression of TDPM genes ", genes[1], ", ", genes[2], " and ", genes[3])))

#Captures last plot

GrabGrob <- function(){
  grid.echo()
  grid.grab()
}

#Creates and stores a correlation plot using corrplot

StoreCorrelationPlot <- function(M, numberCex = 1){
  dev.new(width =  7, height = 7, noRStudioGD = TRUE)
  corrplot(M, method = "color", addCoef.col = "black", tl.col = "black",
           col = wes_palette("Zissou1", 50, type = "continuous"), tl.srt = 45, number.cex = numberCex)
  p <- GrabGrob()
  matrix.colors <- getGrob(p, gPath("square"), grep = TRUE)[["gp"]][["fill"]]
  p <- editGrob(p, gPath("square"), grep = T, gp = gpar(col = NA, fill = NA))
  p <- editGrob(p, gPath("symbols-rect-1"), grep = T, gp = gpar(fill = matrix.colors))
  p <- editGrob(p, gPath("background"), grep = T, gp = gpar(fill = NA))
  p <- grid.arrange(p)
  return (as.ggplot(p))
}

GeneBarPlot <- function(freqTable, yLab, title){
  dev.new(width =  10, height = 10, noRStudioGD = TRUE)
  return(ggplot(freqTable , aes(x = reorder(Gene, Frequency), y = Frequency)) + geom_bar(stat = "identity", fill = "purple") + theme_classic() +
           coord_flip() + geom_text(aes(label = Frequency), hjust = -1, size = 2.5) + xlab("Gene") +
           ylab (yLab) +
           ggtitle(title))
}

MDSPlot <- function(M, nClusters = 6){
  mds.cor <- (1 - M) %>%
    cmdscale() %>%
    as_tibble()
  colnames(mds.cor) <- c("Dim.1", "Dim.2")

  clust <- kmeans(mds.cor, nClusters)$cluster %>%
    as.factor()
  mds.cor <- mds.cor %>%
    ggpubr::mutate(groups = clust)
  # mds.cor <- mds.cor %>%
  #ggpubr::mutate(groups = case_when(
  # clust == 1 ~ "Type 1 stemness-associated clusters",
  # clust == 2 ~ "Type 2 stemness-associated clusters",
  # clust == 3 ~ "Low stemness clusters"))
  p <- ggscatter(mds.cor, x = "Dim.1", y = "Dim.2",
                 label = rownames(M),
                 color = "groups",
                 size = 2,
                 ellipse = T,
                 ellipse.border.remove = T,
                 ellipse.type = "convex",
                 ellipse.alpha = 0.1,
                 font.label = c(10, "bold", "black"),
                 repel = T) + theme(legend.position = "none", legend.title= element_blank())
  return(p)
}

UpsetPlot <- function(upsetInput, upsetName, title, maxDegree = 3, minDegree = 2, minSize = 0, titleHJust = 0,
                      titleVJust = 0){
  return(as.ggplot(upset(
    upsetInput,
    colnames(upsetInput),
    name = upsetName,
    max_degree = maxDegree,
    min_degree = minDegree,
    intersections='all',
    mode = "intersect",
    min_size = minSize,
    max_size = 1000,
    width_ratio = 0.18,
    base_annotations = list(
      'Intersection size' = intersection_size(
        bar_number_threshold = 0,
        mapping=aes(fill='bars_color'),
        mode = "intersect",
        text = list(size = 3.5),
      ) + scale_fill_manual(values=c('bars_color'='navy'), guide='none') +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = titleHJust, vjust = titleVJust, size = 11)) +
        ggtitle(title)
    ),
    set_sizes = upset_set_size(mapping=aes(fill='bars_color')) +
      scale_fill_manual(values=c('bars_color'='navy'), guide = 'none'),
    stripes=c("lightcyan1", "lightblue1"),
    matrix=(
      intersection_matrix(
        geom=geom_point(size=3.5), segment=geom_segment(color ='navy',))
      + scale_color_manual(values=c("white", "navy"), guide = "none")
    ),
  )))
}

CenterTitle <- function(p, textSize = 7)
  return(p + theme(text = element_text(size = textSize), plot.title = element_text(hjust = 0.5)))

NebulosaGrob <- function(seuratObj, features, joint = T, activity = T, newJointTitle = T, palette = "viridis"){
  dev.new(width =  10, height = 10, noRStudioGD = TRUE)
  plots <- lapply(features, function(x) plot_density(seuratObj, x, "data", pal = palette) +
                    theme(text = element_text(size = 7.5), plot.title = element_text(hjust = 0.5)))
  if (joint){
    p <- Plot_Density_Joint_Only(seuratObj, features, viridis_palette = palette) +
      theme(text = element_text(size = 7.5), plot.title = element_text(hjust = 0.5))
    if (newJointTitle)p <- p + ggtitle(str_c("All ",length(features), " genes"))
    plots <- append(plots, list(p))
  }

  if (activity)
    plots <- append(plots, list(plot_density(seuratObj, "activity", "data", pal = palette) +
                                  ggtitle("Activity") + theme(text = element_text(size = 7.5), plot.title = element_text(hjust = 0.5))))
  return(GrobPlot(plots, 3))
}

BarPlotRemovedCells <- function(conditions, nRemovedCells){
  Condition <- unlist(lapply(conditions, function(x) rep(x, 9)))
  Cells <- rep(c("Removed as doublets", "Removed for low novelty score", "Removed as nFeature_RNA outliers",
                 "Removed as nCount_RNA outliers", "Removed for high % of mitochondrial genes",
                 "Removed for high % of ribosomal genes",
                 "Removed for low Shannon diversity", "Removed for low Simpson diversity", "Retained"), length(conditions))
  data <- data.frame(Condition, Cells, nRemovedCells)
  data$Cells <- factor(data$Cells, levels = unique(data$Cells))
  data$Condition <- factor(data$Condition, levels = conditions)
  ggplot(data, aes(fill = Cells, y = nRemovedCells, x = Condition)) +
    geom_bar(position="stack", stat="identity") + theme_classic() +
    scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
                                 '#46f0f0', '#008080', '#e6beff')) +
    ggtitle("Numbers of cells removed at each QC step and retained") + theme(plot.title = element_text(hjust = 0.5)) +
    labs(y = "Number of cells")
}

BarPlotClustersCCSA <- function(clusters, nCCRSA, title = "Numbers of CCRSA markers per cluster"){
  Cluster <- unlist(lapply(clusters, function(x) rep(x, 2)))
  CCRSA <- rep(c("Shared CCRSA markers", "Exclusive CCRSA markers"), length(clusters))
  data <- data.frame(Cluster, CCRSA, nCCRSA)
  data$CCRSA <- factor(data$CCRSA, levels = unique(data$CCRSA))
  data$Cluster <- factor(data$Cluster, levels = clusters)
  ggplot(data, aes(fill = CCRSA, y = nCCRSA, x = Cluster)) +
    geom_bar(position="stack", stat="identity") + theme_classic() +
    scale_fill_manual(values = c('#e6beff', '#e6194b'), name = NULL) +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
    labs(y = "Number of CCRSA markers")
}

BarPlotWilcox <- function(clusters, posNeg){
  Cluster <- unlist(lapply(clusters, function(x) rep(x, 2)))
  nClusters <- rep(c("Number of clusters with lower activity", "- (Number of clusters with higher activity)"), length(clusters))
  data <- data.frame(Cluster, nClusters , posNeg)
  data$nClusters <- factor(data$nClusters , levels = unique(data$nClusters))
  data$Cluster <- factor(data$Cluster, levels = unique(data$Cluster))
  View(data)
  ggplot(data, aes(fill = nClusters, y = posNeg, x = Cluster)) +
    geom_bar(position="stack", stat="identity") + theme_classic() +
    scale_fill_manual(values = wes_palette("Darjeeling1")[c(1, 5)], name = NULL)  +
    ggtitle("Results of pairwise Wilcoxon activity comparisons between clusters") + theme(plot.title = element_text(hjust = 0.5)) +
    labs(y = "Number of clusters")
}

PrepareAlluvialDF <- function(df){
  pvalues <- unique(df$pvalue)
  #Avoiding taking log 0
  if (pvalues[1] == 0)
    df$pvalue <- replace(df$pvalue, df$pvalue == 0, pvalues[2])
  df$Strength <- log(-log(df$pvalue))
  df$pvalue <- c()
  return(df)
}

Alluvialize <- function(df, title, newDev = F, textSize = 3, width = 10, height = 12, titleSize = 12){
  if (newDev)
    dev.new(width = width, height = height, noRStudioGD = TRUE)
  return(ggplot(df,
                aes(y = Strength, axis1 = Cluster, axis2 = Grouping)) +
           geom_alluvium(aes(fill = Grouping), width = 0.5) +
           geom_stratum(aes(fill = Grouping), alpha = 0, width = 0.5) +
           geom_text(stat = "stratum", size = textSize,
                     aes(label = after_stat(stratum))) +
           theme_void() + theme(legend.position = "none") +
           scale_fill_manual(values = distinctColorPalette(length(unique(df$Grouping)))) +
           ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, size = titleSize)))
}

AlluvializeList <- function(dfList, titleList, newDev = F, textSize = 3, titleSize = 12, width = 10, height = 12){
  if (newDev)
    dev.new(width = width, height = height, noRStudioGD = TRUE)
  plots <- mapply(function(x, y) Alluvialize (x, y, textSize = textSize, titleSize = titleSize), dfList, titleList,
                  SIMPLIFY = F)
  return(GrobPlot(plots, length(plots)))
}


IntraclusterCCRSA <- function(df, cond1, cond2, cluster)
  return(ggplot(df, aes(x = avg_log2FC, y = pct.1)) +
           geom_point(color = "purple") + theme_classic() +
           geom_text_repel(aes(label = rownames(df)), size = 3, max.overlaps = Inf) +
           xlab("Average log2 fold change") +
           ylab(str_c("% of ", cond1, " cells expressing gene within ", cluster)) +
           scale_x_continuous(n.breaks = 10) +
           theme(plot.title = element_text(hjust = 0.5)) +
           ggtitle(str_c("Unfiltered CCRSA markers of [", cond1,
                         "] vs. [", cond2, "] within ", cluster)))

PurplePlot <- function(df1, df2, cluster, selection, cutoff = 0){
  shared <- intersect(rownames(df1), rownames(df2))
  df <- data.frame(log1 = df1[shared, ]$avg_log2FC, log2 = df2[shared, ]$avg_log2FC)
  rownames(df) <- shared
  df <- subset(df, log1 > cutoff)
  df <- subset(df, log2 > cutoff)
  View(df)
  return(ggplot(df, aes(x = log1, y = log2)) +
           geom_point(color = "purple") + theme_classic() +
           geom_text_repel(aes(label = rownames(df)), size = 3, max.overlaps = Inf) +
           xlab(str_c("Average log2 fold change - ", cluster)) +
           ylab(str_c("Average log2 fold change - ", selection))  +
           scale_x_continuous(n.breaks = 10) +
           theme(plot.title = element_text(hjust = 0.5)) +
           ggtitle(str_c("Shared markers of the ", cluster, " cluster and the ", selection, "\nselection")))
}

CCRSAPlot <- function(df, width = 12, height = 9){
  dev.new(width = width, height = height, noRStudioGD = TRUE)
  return(ggplot(df, aes(x = node, y = node_degree)) +
           geom_point(color = "purple") + theme_classic() +
           geom_text_repel(aes(label = node), size = 3, max.overlaps = Inf) +
           xlab("CCRSA protein") +
           ylab("Number of interactions in the CCRSA protein-protein interaction network")  +
           ggtitle("Number of interactions in the CCRSA protein-protein interaction network") +
           theme(plot.title = element_text(hjust = 0.5),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank()))
}

