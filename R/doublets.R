#' @export

PredictDoublets <- function(seuratObj, start = 1, stop = 100)
{
  #Adds an user-defined number of scDblFinder runs to the Seurat metadata
  DefaultAssay(seuratObj) <- "RNA"
  sce_seurat <- as.SingleCellExperiment(seuratObj)
  for (i in start:stop){
    stri <- toString(i)
    message(str_c("Predicting doublets: run ", stri, "."))
    seuratObj[[str_c("unit.class", stri)]] <- scDblFinder(sce_seurat, samples = "orig.ident")$scDblFinder.class
    gc()
  }
  return(seuratObj)
}

AddDoubletRuns <- function(seuratObj, datasetName, nRuns = 100, nSteps = 4){
  #Also saves the metadata after adding the runs
  #Using 4 steps by default to avoid memory issues
  stepSize <- as.integer(nRuns/nSteps)
  for (i in 1: (nSteps - 1))seuratObj <- PredictDoublets(seuratObj, (i - 1) * stepSize + 1, i * stepSize)
  seuratObj <- PredictDoublets(seuratObj, (nSteps - 1) * stepSize + 1, nRuns)
  write.csv(seuratObj@meta.data, str_c(datasetName, " ", nRuns, " doublet runs metadata.csv"))
  return(seuratObj)
}

AddDoubletListings <- function(seuratObj, nRuns = 100){
  #Counts how mnay time each cell was predicted to be a doublet and adds this information to the Seurat metadata
  isDoublet <- as.data.frame(lapply(1:nRuns,
                                    function(x)str_count(seuratObj@meta.data[[str_c("unit.class", x)]], "doublet")))
  seuratObj$doublet.listings <- rowSums(isDoublet)
  return(seuratObj)
}

FindMeanDoublets <- function(seuratObj, nRuns = 100)
  return(sum(seuratObj$doublet.listings)/nRuns)

FindDoubletCounts <- function(seuratObj, nRuns = 100)
  #Starting from 2 in order to exclude cells never predicted to be doublets from the result
  return(dplyr::count(seuratObj@meta.data, doublet.listings)$n[2:(nRuns + 1)])

FindReverseCumSum <- function(doubletCounts)
  return(rev(cumsum(rev(doubletCounts))))

FindDoubletCutoff <- function(meanDoublets, reverseCumSum){
  #A doublet.listings column must have been added to the Seurat object
  message(str_c("The average number of doublets predicted per run was ",  meanDoublets))
  cutoff <- which.min(abs(reverseCumSum - meanDoublets))
  message(str_c("Cell units predicted to be doublets in at least ",  cutoff, " runs will be regarded as doublets."))
  message(str_c("Taking this cutoff implies ", reverseCumSum[cutoff], " predicted doublets."))
  return(cutoff)
}

AcceptedDoubletsPIDF <- function(doubletCounts, cutoff, reverseCumSum, meanDoublets, nRuns = 100){
  df <- data.frame(Listings = 1:nRuns, Count = doubletCounts)
  df <- mutate(df, Outcome = if_else(Listings >= cutoff, "Accepted doublets", "Rejected doublets"))
  df$arcs <- reverseCumSum - meanDoublets
  return(df)
}

#originally: width = 8.1, height = 9
PlotAcceptedDoublets <- function(doubletCounts, cutoff, reverseCumSum, meanDoublets, nRuns = 100,
                                 width = 8.1, height = 9, textSize = 7){
  df <- AcceptedDoubletsPIDF(doubletCounts, cutoff, reverseCumSum, meanDoublets)
  dev.new(width =  width, height = height, noRStudioGD = TRUE)
  plots <- list(ggplot(df, aes(x = Listings, y = arcs)) +
                  geom_point(stat = "identity", color = c(rep("red", cutoff - 1),
                                                          "purple", rep("blue", nRuns - cutoff))) +
                  ggtitle("Identification of doublet prediction significance cutoff") +
                  xlab("Starting index of reverse cumulative sum") +
                  ylab("Difference between reverse cumulative sum\nand mean number of doublets") +
                  geom_vline(aes(xintercept = cutoff), color="purple", size = 0.5) +
                  geom_hline(aes(yintercept = 0), color="pink", size = 0.5) +
                  theme_classic() +
                  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = textSize + 2)) +
                  geom_text(data = data.frame(cutoff), aes(x = cutoff, label = str_c("x = ", cutoff)),
                            y = 150, angle = -10, size = 3),
                ggplot(df, aes(x = Listings, y = Count, fill = Outcome)) +
                  geom_bar(stat = "identity") +
                  scale_fill_manual(values = c("navy", "aliceblue")) +
                  ggtitle("Predicted doublets after 100 scDblFinder runs") +
                  xlab("Percentage of runs") +
                  ylab("Number of cells predicted to be doublets") +
                  theme_classic() +
                  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = textSize + 2))
  )
  return(GrobPlot(plots, 2, topMargin = 5, bottomMargin = 5))
}

AddConsensusDoublets <- function(seuratObj, cutoff){
  seuratObj@meta.data <- seuratObj@meta.data %>%
    mutate(unit.class = case_when(
      doublet.listings < cutoff ~ "singlet",
      doublet.listings >= cutoff ~ "doublet"
    ))
  #Factoring is necessary to assure correct order ((singlet, doublet))
  seuratObj$unit.class <- factor(x = seuratObj$unit.class, levels = c("singlet", "doublet"))
  return(seuratObj)
}

Jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

RunsJaccards <- function(runs, nRuns = 100)
  return(lapply(1:nRuns, function(x)
    sort(sapply(1:nRuns, function(y) Jaccard(runs[[x]], runs[[y]])))[1:(nRuns - 1)]))

DoubletsJaccardsPIDF <- function(seuratObj, nRuns = 100){
  consensus <- which(seuratObj[["unit.class"]] == "doublet")
  runs <- lapply(1:nRuns, function(x) which(seuratObj[[str_c("unit.class", x)]] == "doublet"))
  runsJaccards <- RunsJaccards(runs, nRuns)
  consensusJaccards <- sapply(1:100, function(x) Jaccard(consensus, runs[[x]]))
  maxJaccards <- sapply(runsJaccards, function(x) max(x))
  meanJaccards <- sapply(runsJaccards, function(x) mean(x))
  minJaccards <- sapply(runsJaccards, function(x) min(x))
  df <- data.frame(Index = 1:nRuns,
                   Consensus = consensusJaccards,
                   Maximum = maxJaccards,
                   Mean = meanJaccards,
                   Minimum = minJaccards)
  return(df)
}

PlotDoubletsJaccards <- function(seuratObj, nRuns = 100, width = 8.1, height = 9, textSize = 7){
  df <- DoubletsJaccardsPIDF(seuratObj, nRuns)
  dev.new(width =  width, height = height, noRStudioGD = TRUE)
  ggplot(df, aes(x = Index, y = Consensus)) +
    geom_point(color = "purple", shape = 2) +
    geom_point(aes(y = Maximum), color = "navy", shape = 0) +
    geom_point(aes(y = Mean), color = "royalblue", shape = 1) +
    geom_point(aes(y = Minimum), color = "slategray1", shape = 6) +
    geom_hline(yintercept = mean(df$Consensus), color = "purple", size = 1, linetype = "dashed") +
    geom_hline(yintercept = mean(df$Maximum), color = "navy", size = 1, linetype = "dashed") +
    geom_hline(yintercept = mean(df$Mean), color = "royalblue", size = 1, linetype = "dashed") +
    geom_hline(yintercept = mean(df$Minimum), color = "slategray1", size = 1, linetype = "dashed") +
    ylim(0, 1) +
    ggtitle("Jaccard similarity scores of pairs of scDblFinder runs and the consensus prediction") +
    xlab("scDblFinder run index") +
    ylab("Jaccard similarity score") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = textSize + 2))
}

StripDoubletRuns <- function(seuratObj, datasetName, nRuns = 100){
  #Also saves the metadata after stripping the runs
  for (i in 1:nRuns) seuratObj@meta.data[[str_c("unit.class", i)]] <- c()
  write.csv(seuratObj@meta.data, str_c(datasetName, " stripped metadata.csv"))
  return(seuratObj)
}

ImportDoublets <- function(seuratObj, datasetName){
  seuratObj$unit.class <- read.csv(str_c(datasetName, " stripped metadata.csv"), row.names = 1)$unit.class
  seuratObj$unit.class <- factor(x = seuratObj$unit.class, levels = c("singlet", "doublet"))
  return(seuratObj)
}

AddDoubletInformation <- function(seuratObj, datasetName, doFindDoublets, doPlotDoublets, nRuns = 100){
  if (doFindDoublets)seuratObj <- AddDoubletRuns(seuratObj, datasetName)
  if (doPlotDoublets){
    seuratObj <- AddDoubletListings(seuratObj, nRuns)
    meanDoublets <- FindMeanDoublets(seuratObj, nRuns)
    doubletCounts <- FindDoubletCounts(seuratObj, nRuns)
    reverseCumSum <- FindReverseCumSum(doubletCounts)
    doubletCutoff <- FindDoubletCutoff(meanDoublets, reverseCumSum)
    print(PlotAcceptedDoublets(doubletCounts, doubletCutoff, reverseCumSum, meanDoublets, nRuns))
    seuratObj <- AddConsensusDoublets(seuratObj, doubletCutoff)
    print(PlotDoubletsJaccards(seuratObj, nRuns))
    seuratObj <- StripDoubletRuns(seuratObj, datasetName, nRuns)
  }else seuratObj <- ImportDoublets(seuratObj, datasetName)
  return(seuratObj)
}

RestoreRuns <- function(seuratObj, datasetName, nRuns = 100){
  handy <- read.csv(str_c(datasetName, " ", nRuns, " doublet runs metadata.csv"))
  for (i in 1:nRuns)
    seuratObj@meta.data[[str_c("unit.class", i)]] <- handy[[str_c("unit.class", i)]]
  return(seuratObj)
}

