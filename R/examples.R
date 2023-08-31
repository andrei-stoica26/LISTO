#Commented code. Various bits exemplifying how functions from the package are used.

# seuratObj <- CreateSeuratObject(counts = readRDS(files[1]), project = conditions[1], min.cells = 10)
#
#
# seuratObj <- MergeSeurats(files = str_c("PatientData/", list.files("PatientData/")),
#                           conditions = conditions,
#                           ids = c("1D", "2I", "3DG", "4IG"))

#

# doFindDoublets <- F
# seuratObj <- MergeSeurats(files, conditions, ids)
# if (doFindDoublets)
#   seuratObj <- AddDoubletInformation(seuratObj, "Patient", T, T) else{
#     seuratObj <- RestoreRuns(seuratObj, "Patient")
#     seuratObj <- AddDoubletInformation(seuratObj, "Patient", F, T)
#   }
#
# seuratObj <- AddFilteringCriteria(seuratObj, T)

#
# seuratObj <- QualityControl(seuratObj,
#                             doRemoveDoublets = T,
#                             noveltyBounds = lapply(1:length(conditions), function(x) allBounds[[x]][[2]]),
#                             featureBounds = lapply(1:length(conditions), function(x) allBounds[[x]][[3]]),
#                             countBounds = lapply(1:length(conditions), function(x) allBounds[[x]][[4]]),
#                             mtBounds = lapply(1:length(conditions), function(x) allBounds[[x]][[5]]),
#                             riboBounds = lapply(1:length(conditions), function(x) allBounds[[x]][[6]]),
#                             shannonBounds = lapply(1:length(conditions), function(x) allBounds[[x]][[7]]),
#                             simpsonBounds = lapply(1:length(conditions), function(x) allBounds[[x]][[8]]))
# rareGenes <- FindRareGenes(seuratObj)
# seuratObj <- RemoveRareGenes(seuratObj, rareGenes)
# seuratObj <- IndividualSCT(seuratObj)
# seuratObj <- RunPCA(seuratObj)
# seuratObj <- RunHarmony(seuratObj, assay = "SCT", "orig.ident")
# pct <- seuratObj[["harmony"]]@stdev / sum(seuratObj[["harmony"]]@stdev) * 100
# nUMAPDims <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# seuratObj <- RunUMAP(seuratObj, reduction = "harmony", dims = 1:nUMAPDims)
# seuratObj <- FindNeighbors(seuratObj, reduction = "umap", dims = 1:2)
# seuratObj <- FindClusters(seuratObj, resolution = 0.03)

# M <- GetGeneCorrelationMatrix(seuratObj, geneSets$trad)
# View(M)
# StoreCorrelationPlot(M, 0.8) + ggtitle("TDPM genes - Spearman correlation plot of expression") +
#   theme(plot.title = element_text(hjust = 0.5))
# View(M)
#
# geneSets <- AddGeneSets(seuratObj, tradMarkersPatient, LitMarkers, sideLogMarkers)
#

# x <- LitVLogV(geneSets$Lit, upMarkersPatient, levels(seuratObj), length(rownames(seuratObj)))
# x[[14]]
#
#
# prepbarCCRSA <- CountCCRSA(upNamesPatient, geneSets)
# BarPlotClustersCCSA(0:(length(upNamesPatient) - 1), prepbarCCRSA)
#
# x <- LogLogV(geneSets$side, upMarkersPatient, levels(seuratObj), length(rownames(seuratObj)))
# index <- 7
# m <- upMarkersPatient[[index]][intersect(rownames(geneSets$side), upNamesPatient[[index]]), ]
# m <- m[order(m$p_val_adj),]
# MessageVector(rownames(m))
# print(m$p_val_adj)
#
#
# df <- TopActivityGroupings(seuratObj, seurat_clusters, levels(seuratObj), F)
#
# PrintOverlap(df)
# levels(seuratObj$orig.ident)
#
#
# df2 <- ActivityWilcoxCluster(seuratObj)
# sc <- ScoreActivityWilcox(seuratObj, df2)
# posNeg <- unlist(mapply(function(x, y) c(x, -y), sc$Wins, sc$Losses, SIMPLIFY = F))
# BarPlotWilcox(gsub("Cluster ", "", sc$Cluster), posNeg)

# curves <- readRDS("CurvesPatient.rds")
# counts <- seuratObj@assays$SCT@counts
# gamPatient <- readRDS("PatientTradeseqGam.rds")
# PlotLineages(seuratObj, curves, dimred)
# seuratObj <- AddLineages(seuratObj, curves)
# seuratObj <- AddCurveweights(seuratObj, curves)
# plotGeneCount(curve = curves, counts = counts,
#               clusters = apply(slingClusterLabels(curves), 1, which.max),
#               models = gamPatient) + ggtitle("Differentiation lineages") + theme(plot.title = element_text(hjust = 0.5))

# LineagePlot(seuratObj, dimred, curves, "Lineage 1")
# LineagePlot(seuratObj, dimred, curves, "Lineage 2")
# LineagePlot(seuratObj, dimred, curves, "Lineage 3")
#
# enrichSets <- c("CellMarker_Augmented_2021",
#                 "PanglaoDB_Augmented_2021")
# enriched <- lapply(upNamesPatient, function(x) rba_enrichr(gene_list = x, gene_set_library = enrichSets))
# cellMarker <- EnrichrVToER(enriched, "CellMarker_Augmented_2021")
# panglao <- EnrichrVToER(enriched, "PanglaoDB_Augmented_2021")
# ConceptCouple(cellMarker, 4, 7, 6, "CellMarker_Augmented_2021")
# ConceptSolo(cellMarker, 8, 7, "CellMarker_Augmented_2021")
# ConceptSolo(panglao, 4, 1, "PanglaoDB_Augmented_2021")
#
# View(ClusterGO)
# clusterEnrichments <- ClusterGO(seuratObj, upMarkersPatient)
# clusterDescriptions <- EnrichmentDescriptions(clusterEnrichments)
# saveRDS(clusterEnrichments, "PatientClusterEnrichments.rds")
# saveRDS(clusterDescriptions, "PatientClusterDescriptions.rds")
#
# options(ggrepel.max.overlaps = Inf)
# ConceptCouple(clusterEnrichments, 4, 7, 3, "GO")
# ConceptCouple(clusterEnrichments, 1, 5, 3, "GO")
# ConceptCouple(clusterEnrichments, 2, 3, 2, "GO")
# ConceptCouple(clusterEnrichments, 6, 8, 2, "GO")
# ConceptSolo(clusterEnrichments, 0, 3, "GO")
#
# seuratObj <- AnnotateSeurat(seuratObj, "cell.type",
#                             c("Actin\norganization", "Hypoxia-like", "Stress\nresponse",
#                               "Protein\ndegradation", "Top\nstemness", "Cellular\nrespiration",
#                               "Wnt\nsignalling", "p53\nsignalling", "IF\nresponse"))
# DimPlot(seuratObj, group.by = "cell.type", label = T) + ggtitle("The functional annotation of the clusters") +
#   theme(plot.title = element_text(hjust = 0.5))
#
# length(rownames(seuratObj))
#
# gamPatient <- readRDS("PatientTradeseqGam.rds")
# counts <- seuratObj@assays$SCT@counts
# plotSmoothers(gamPatient, counts, "BIRC5", lwd = 0.8)
# SaveEarlies(gamPatient, 10, "Patient")
# earlies <- readRDS("EarliesPatient.rds")
# SaveStartends(gamPatient, "Patient")

# cellchat <- readRDS("PatientCellchat.rds")
# cellchat@netP$pathways
#
# NetVisualAgg(cellchat, cellchat@netP$pathways[51:58])
#
# NetVisualAgg(cellchat, c("MK","EGF", "MIF", "THY1", "PTN", "EDA", "CLDN", "PECAM1"))
# NetVisualAgg(cellchat, c("DESMOSOME", "PDGF", "PTPRM", "PROS", "SPP1", "SEMA6", "IFN-I", "JAM"))
# p1 <- selectK(cellchat, pattern = "outgoing")
# p2 <- selectK(cellchat, pattern = "incoming")
#
# p2
#
# cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 4)
# cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 6)
#
# netAnalysis_river(cellchat, pattern = "outgoing")
# netAnalysis_river(cellchat, pattern = "incoming")
# dev.new(width =  12, height = 12, noRStudioGD = TRUE)
# GrobPlot(list(netAnalysis_dot(cellchat, pattern = "outgoing", font.size = 6),
#               netAnalysis_dot(cellchat, pattern = "incoming", font.size = 6)), 2)

# ##################################Epigenetics
# allEpi <- "histone|chromatin|epigenetic|protein-DNA|packaging|nucleosome|
# DNA methylation|imprinting|silencing|ncRNA-mediated|conformation|geometric|miRNA|inactivation of X chromosome"
# pathways <- FindPathways(seuratObj)
# epiPathways <- FindEpigeneticPathways(pathways, allEpi)
# clusterDescriptions <- readRDS("PatientClusterDescriptions.rds")
# clEpiPathways <- FindGroupingsEpiPathways(clusterDescriptions, epiPathways)
# Reduce(union, clEpiPathways)
# DimPlot(seuratObj, group.by = "cell.type", label = T)
# types <-  gsub("\n", " ", c("Actin\norganization", "Hypoxia-like", "Stress\nresponse",
#                             "Protein\ndegradation", "Top\nstemness", "Cellular\nrespiration",
#                             "Wnt\nsignalling", "p53\nsignalling", "IF\nresponse"))
# CountGroupingEpiPathways(clEpiPathways, types)
# PrintGroupingEpiPathways(clEpiPathways, types)
# clusterEnrichments <- readRDS("PatientClusterEnrichments.rds")
# EpigeneticSolo(clusterEnrichments, 5, epiPathways, types, "cluster")
# LitLogV(epiPathways, clusterDescriptions, types, length(pathways), F)
#
# ################################Pseudobulk
# DimPlot(seuratObj, label = T)
# activityStats4 <- ActivityStats(RemoveNullGenes(subset(seuratObj, seurat_clusters == 4)))
# saveRDS(activityStats4, "PatientActivityStats4.rds")
# activityStats7 <- ActivityStats(RemoveNullGenes(subset(seuratObj, seurat_clusters == 7)))
# saveRDS(activityStats7, "PatientActivityStats7.rds")
# slingshotStats4 <- SlingshotStats(RemoveNullGenes(subset(seuratObj, seurat_clusters == 4)))
# saveRDS(slingshotStats4, "PatientSlingshotStats4.rds")
# slingshotStats7 <- SlingshotStats(RemoveNullGenes(subset(seuratObj, seurat_clusters == 7)))
# saveRDS(slingshotStats7, "PatientSlingshotStats7.rds")

# View(subset(activityStats4, NCells > 9))
# View(subset(slingshotStats4, NCells > 9))
# View(subset(activityStats7, NCells > 9))
# View(subset(slingshotStats7, NCells > 9))
# View(activityStats7)
# dim(seuratObj)
#
# subset(seuratObj, seurat_clusters == 4)$activity
# View(activityStats4)
#
# sigTrad <- c("EZH2", "NES", "ABCG2", "ABCB1")
# SignatureRepresentation(seuratObj, sigTrad)
# sigTrad2 <- c("EZH2", "NES", "THY1")
# SignatureRepresentation(seuratObj, sigTrad2)
# Plot_Density_Joint_Only(seuratObj, sigTrad, "viridis") + theme(plot.title  = element_text(hjust = 0.5))
# Plot_Density_Joint_Only(seuratObj, sigTrad2, "viridis") + theme(plot.title  = element_text(hjust = 0.5))
# selectionMarkers <- readRDS("PatientSelectionMarkers.rds")
# invisible(lapply(LitMarkersOverlap(seuratObj, geneSets, selectionMarkers, RownamesCMLPatient), PrintOverlap))
# PrintOverlap(SideMarkersOverlap(seuratObj, geneSets, selectionMarkers, RownamesCMLPatient))
#
# PrintOverlap(TopActivityGroupings(seuratObj, orig.ident, levels(seuratObj$orig.ident), F))
# PrintOverlap(TopActivityGroupings(seuratObj, orig.ident, levels(seuratObj$orig.ident), T))
# PrintOverlap(ActivityWilcoxOrig(seuratObj))
#
# ##################################Intracluster
# PrintOCC(CheckCondDist(seuratObj, F), types)
# PrintOCC(CheckCondDist(seuratObj, T), types)
# PrintOCC(CheckCondDist(subset(seuratObj, orig.ident %in% levels(seuratObj$orig.ident)[1:2]), F), types)
#
# View(CheckCondDist(seuratObj, F))
#
# seurats <- OrderedSplit(seuratObj)
# SignatureRepresentation(seurats[[2]], sigTrad2)
# types
# ICSMarkers <- FindICSMarkers(seurats, RownamesCMLPatient, firstsPatient, secondsPatient)
# saveRDS(ICSMarkers, "ICSMarkersPatient.rds")
# ICSMarkers <- readRDS("ICSMarkersPatient.rds")
#
# invisible(lapply(1:14, function(index)
#   PrintOverlap(
#     LitSelectionsInClusters(geneSets$Lit, ICSMarkers, RownamesCMLPatient, index,  types, length(rownames(seuratObj))))
# ))
# PrintOverlap(
#   LogSelectionsInClusters(rownames(geneSets$side), ICSMarkers, RownamesCMLPatient,  types, length(rownames(seuratObj))))
#
# wilcoxClusters <- ActivityWilcoxOrigWithinClusters(seurats, seuratObj, types)
# alluvialDF <- PrepareAlluvialDF(wilcoxClusters)
# dev.new(width =  12, height = 12, noRStudioGD = TRUE)
# Alluvialize(alluvialDF, "Intracluster activity Wilcoxon pairwise comparisons")
#
# ActivityCohen(seurats[[2]], "DMSO", "I-BRD9")


#
# ###########################Conditions and clusters
# CCMarkers <- CCOverlap(seuratObj, upMarkersPatient, selectionMarkers, RownamesCMLPatient, types, length(rownames(seuratObj)),
#                        isLog = T)


# TwoLogsJaccard(upMarkersPatient[[5]], selectionMarkers[[8]], length(rownames(seuratObj)))
#
# alluvialCC <- PrepareAlluvialDF(CCMarkers[1:20, ])
# Alluvialize(alluvialCC, "Top 20 overlaps between cluster markers and selection markers")

# selectionEnrichments <- SelectionGO(seuratObj, selectionMarkers)
# selectionDescriptions <- EnrichmentDescriptions(selectionEnrichments)
#
#
# nPathways <- length(rownames(GenesER(rownames(seuratObj))@result))
#
# CCPathwaysMarkers <- CCPathwaysOverlap(seuratObj, clusterDescriptions, selectionDescriptions,
#                                        RownamesCMLPatient, types, nPathways)
#
#
# alluvialP <- PrepareAlluvialDF(CCPathwaysMarkers[1:20, ])
# Alluvialize(alluvialP, "Top 20 overlaps between GO terms enriched for clusters and GO terms enriched for selections")
#
#
# x <- AllLitCCOverlap(geneSets$Lit, upMarkersPatient, selectionMarkers,
#                      RownamesCMLPatient, types, length(rownames(seuratObj)))
# saveRDS(x, "ThreeWayLitPatient.rds")
#
# PrintOverlapCC(LogCCOverlap(geneSets$side, upMarkersPatient, selectionMarkers, RownamesCMLPatient,
#                             types, length(rownames(seuratObj))))
#
#
# rm(alluvial3)
# alluvial3 <- PrepareAlluvialDF(x[[14]])
# AlluvializePair(PrepareAlluvialDF(x[[2]]),
#                 PrepareAlluvialDF(x[[7]]),
#                 "Overlaps between cluster markers, selection markers, and the pancreas CCRSA set",
#                 "Overlaps between cluster markers, selection markers, and the glioma CCRSA set")



# allEpi <- "histone|chromatin|epigenetic|protein-DNA|packaging|nucleosome|
# DNA methylation|imprinting|silencing|ncRNA-mediated|conformation|geometric|miRNA|inactivation of X chromosome"
# pathways <- FindPathways(seuratObj)
# selectionDescriptions <- readRDS("PatientSelectionDescriptions.rds")
# length(epiPathways)
# epiPathways <- FindEpigeneticPathways(pathways, allEpi)
# selEpiPathways <- FindGroupingsEpiPathways(selectionDescriptions, epiPathways)
# length(Reduce(union, selEpiPathways))
# CountGroupingEpiPathways(selEpiPathways, RownamesCMLPatient)
# PrintOverlap(LitLogV(epiPathways, clusterDescriptions, types, length(pathways), F))
#
# RownamesCMLPatient
# EpigeneticSolo(selectionEnrichments, 3, epiPathways, RownamesCMLPatient, "selection")
# EpigeneticSolo(selectionEnrichments, 2, epiPathways, RownamesCMLPatient, "selection")
# #Section: Further BRD9
# selectionNames <- lapply(selectionMarkers, rownames)

# enrichSets <- c("ESCAPE", "KEA_2015", "Gene_Perturbations_from_GEO_up", "Jensen_DISEASES", "ARCHS4_TFs_Coexp",
#                 "Enrichr_Submissions_TF-Gene_Coocurrence", "CellMarker_Augmented_2021",
#                 "PanglaoDB_Augmented_2021")
# enriched <- lapply(upNames, function(x) rba_enrichr(gene_list = x, gene_set_library = enrichSets))
#
# allEnrichSets <- names(enriched[[7]])
# which (allEnrichSets %in% "PanglaoDB_Augmented_2021")
#
# ViewEnriched <- function(index, cluster){
#   print(allEnrichSets[index])
#   View(enriched[[cluster + 1]][[allEnrichSets[index]]])
# }
#
# jensen <- EnrichrVToER(enriched, "Jensen_DISEASES")
#
# #Section 6: Trajectory
# gamTradeseq <- readRDS("A13ATradeseqGam.rds")
# counts <- readRDS("CountsA13A.rds")
# curves <- readRDS("CurvesA13A.rds")
# earlies <- readRDS("EarliesA13A.rds")
#
# plotGeneCount(curve = curves, counts = counts,
#               clusters = apply(slingClusterLabels(curves), 1, which.max),
#               models = gamTradeseq) + ggtitle("Differentiation lineages") +
#   theme(plot.title = element_text(hjust = 0.5))
#
# earlyStage <- c("FTL","AKR1C1","AKR1C2", "NQO1", "TXNRD1", "TALDO1")
# middleStage <- c("AKR1B10", "ASPH", "PALS2", "EZR", "TNFRSF11B", "GDF15")
# lateStage <- c("CEACAM6", "REG4", "ANXA1", "TSPAN8", "ANKRD1", "GPRC5A")
#
# SmoothersGrob(gamTradeseq, counts, earlyStage, "Genes marking the separation of Lineage 3")
# SmoothersGrob(gamTradeseq, counts, middleStage, "Genes marking the middle course of Lineage 3")
# SmoothersGrob(gamTradeseq, counts, lateStage,"Genes marking the late stage of Lineage 3")
#
# ueg <- FindUniqueEarliesGenes(earlies, upNames, 3, fc = 0.1, wald = 20)
#
# separation3 <- earlyDETest(gamTradeseq, knots = c(4, 6))
# sepRestricted <- subset(separation3[intersect(rownames(separation3), upNames[[3]]),], pvalue < 0.05 &
#                           fcMedian > 0.1 & waldStat > 20)
# sepRestricted <- sepRestricted[order(sepRestricted$waldStat, decreasing = T),]
# separationGenes <- rownames(sepRestricted)
# paste(separationGenes, collapse =", ")
# MessageVector(separationGenes)
# sepGO <- GenesER(rownames(sepRestricted))
# sepGO@result <- subset(sepGO@result, p.adjust < 0.05)
# sepWP <- GenesER(rownames(sepRestricted), "enrichWP")
# sepWP@result <- subset(sepWP@result, p.adjust < 0.05)
# dev.new(width =  12, height = 12, noRStudioGD = TRUE)
# GrobPlot(list(StoreCnetplot(sepGO, 20), StoreCnetplot(sepWP, 7)), 2, leftMargin = 5) +
#   ggtitle(str_c("GO (A) and WP (B) concept network plot of the early drivers of the separation of Lineage 3")) +
#   theme(plot.title = element_text(hjust = 0.5))
#
# late3 <- earlyDETest(gamTradeseq, knots = c(9, 11))
# latRestricted <- subset(late3[intersect(rownames(late3), upNames[[3]]),], pvalue < 0.05 &
#                           fcMedian > 0.1 & waldStat > 20)
# MessageVector(rownames(latRestricted))
# latWP <- GenesER(rownames(latRestricted), "enrichWP")
# latWP@result <- subset(latWP@result, p.adjust < 0.05)
# dev.new(width =  6, height = 6, noRStudioGD = TRUE)
# StoreCnetplot(latWP) + ggtitle(str_c("WP concept network plot of markers of the late evolution of Lineage 3")) +
#   theme(plot.title = element_text(hjust = 0.5, size = 10))
#
# RedYellowPlot(curves, counts, c("ANKRD1", "VSIG1", "GPRC5A"),
#               "Feature plots of SP genes on the direction of nascent RNA from Transitional SP")
#
# seuratObj <- RunMonocle(seuratObj, upMarkers[[7]])
#
# activityStats <- readRDS("ActivityStatsA13A.rds")
# slingshotStats <- readRDS("SlingshotStatsA13A.rds")
# monocleStats <- readRDS("MonocleStatsA13A.rds")
#
# View(activityStats)
#
# StemnessExtremaGrob(activityStats, slingshotStats, monocleStats)
#
# ActivityLPG(activityStats)
# SlingshotLPG(slingshotStats)
# MonocleLPG(monocleStats)
#
# lpGenes <- FindLocalPeakGenes(activityStats, slingshotStats, monocleStats)
#
# lpFreqTable <- LocalPeaksFrequencyTable(lpGenes)
#
# paste(setdiff(lpFreqTable$Gene, geneSets$Lit[[14]]), collapse = ", ")
# MessageVector(setdiff(lpFreqTable$Gene, c("DEPDC1", geneSets$Lit[[14]])))
# MessageVector(intersect(lpFreqTable$Gene, geneSets$Lit[[14]]))
#
# View(GenesER(c(setdiff(lpFreqTable$Gene, c(geneSets$Lit[[14]], "DEPDC1"))), "enrichGO")@result)
#
# View(GenesER(c(setdiff(lpFreqTable$Gene, geneSets$Lit[[14]])))@result)
#
# GeneBarPlot(freqTable, "Number of occurrences among the 86 genes linked with median and mean ORIGINS activity,
#         Slingshot pseudotime and Monocle 3 pseudotime",
#             "Bar plot of the genes with at least 3 stemness assay median and mean associations")
# LPUpsetGrob(lpGenes)
# PrintOverlap(MedianLPGenesOverlap(seuratObj, lpGenes))
# PrintOverlap(MeanLPGenesOverlap(seuratObj, lpGenes))
#
# PrintOverlap2(LitLitV(geneSets$Lit[[14]], lpGenes, GetPeakSetLabels(), nrow(seuratObj)))
#
# NebulosaGrob(seuratObj, c("HSP90AB1", "PTMA", "STMN1", "PARP1", "NUCKS1", "HSP90AA1", "HMGN2"))
# NebulosaGrob(seuratObj, c("HMGB2", "UBB", "SCG2", "RPS6", "HNRNPA2B1","H2AZ1",  "DEPDC1"))
#
# #Section 7: Pseudobulk
# o <- readRDS("ActivityConditionStatsA13A.rds")
# p <- readRDS("SlingshotConditionStatsA13A.rds")
# q <- readRDS("MonocleConditionStatsA13A.rds")
#
#
# sigTrad <- c("AGR2", "ALDH1A1", "REG4", "TSPAN8")
# SignatureRepresentation(seuratObj, sigTrad)
# SignatureRepresentation(subset(seuratObj, orig.ident != "SB-431542"), sigTrad)
# SignatureRepresentation(subset(seuratObj, seurat_clusters == 2), sigTrad)
# SignatureRepresentation(subset(seuratObj, orig.ident != "SB-431542" & seurat_clusters == 2), sigTrad)
# Plot_Density_Joint_Only(seuratObj, sigTrad, "viridis") + theme(plot.title  = element_text(hjust = 0.5))
# VlnPlot(seuratObj, "activity", group.by = "orig.ident")
#
#
#
# PrintAllSGSelections(geneSets$trad, selectionMarkers, RownamesCMLA13A)
#
# invisible(lapply(LitMarkersOverlap(seuratObj, geneSets, selectionMarkers, RownamesCMLA13A), PrintOverlap))
# PrintOverlap(SideMarkersOverlap(seuratObj, geneSets, selectionMarkers, RownamesCMLA13A))
#
# PrintOverlap(TopActivityGroupings(seuratObj, orig.ident, levels(seuratObj$orig.ident)))
# PrintOverlap(TopActivityGroupings(seuratObj, orig.ident, levels(seuratObj$orig.ident), T))
#
# acts <- rownames(GetActivityCorrelations(seuratObj))[1:50]
# PrintOverlap(LitLogV(acts, selectionMarkers, RownamesCMLA13A, length(rownames(seuratObj))))
# invisible(lapply(1:length(selectionMarkers), function(i){
#   print(RownamesCMLA13A[i])
#   print(paste(intersect(rownames(selectionMarkers[[i]]), acts), collapse = ", "))
# }))
# PrintOverlap(ActivityWilcoxOrig(seuratObj))
# PrintOverlap(SlingshotWilcoxOrig(seuratObj))
#
# ActivityHist(seuratObj)
#
# TopSlingshotGroupings(seuratObj, orig.ident, levels(seuratObj$orig.ident), lowerTail = F)
# TopMonocleGroupings(seuratObj, orig.ident, levels(seuratObj$orig.ident), lowerTail = F)
# TopActivityGroupings(seuratObj, orig.ident, levels(seuratObj$orig.ident), lowerTail = F)
#
# #Section 8: Intracluster
# PrintOCC(CheckCondDist(seuratObj, F), types)
# PrintOCC(CheckCondDist(seuratObj, T), types)
# PrintOCC(CheckCondDist(subset(seuratObj, orig.ident %in% levels(seuratObj$orig.ident)[1:2]), F), types)
#
# CCMarkers <- readRDS("CCMarkersA13A.rds")
# invisible(lapply(1:length(CCMarkers), function(i){
#   message(types[i])
#   PrintAllSGSelections(geneSets$trad, CCMarkers[[i]], RownamesCMLA13A)
# }))
#
# invisible(lapply(1:14, function(index)
#   PrintOverlap(
#     LitSelectionsInClusters(geneSets$Lit, CCMarkers, RownamesCMLA13A, index,  types, length(rownames(seuratObj))))
# ))
#
# PrintOverlap(
#   LogSelectionsInClusters(rownames(geneSets$side), CCMarkers, RownamesCMLA13A,  types, length(rownames(seuratObj))))
#
# seurats <- OrderedSplit(seuratObj)
# PrintOverlapCC(ActivityWilcoxOrigWithinClusters(seurats, seuratObj, types))
# PrintOverlapCC(SlingshotWilcoxOrigWithinClusters(seurats, seuratObj, types))
# PrintOverlapCC(MonocleWilcoxOrigWithinClusters(seurats, seuratObj, types))
#
# hist(subset(seuratObj, seurat_clusters == 3)$monocle)
#
# hist(subset(seuratObj, seurat_clusters == 0)$slingshot)
#
# DimPlot(seuratObj, label = T)
# VlnPlot(seurats[[7]], "monocle", group.by = "orig.ident")
#
# #Section 9: Clusters and conditions
# selectionMarkers <- readRDS("A13ASelectionMarkers.rds")
# PrintOverlapCC(CCOverlap(seuratObj, upMarkers, selectionMarkers, RownamesCMLA13A, types, length(rownames(seuratObj))))
#
# upDescriptions <- readRDS("UpDescriptionsA13A.rds")
# selectionDescriptions <- readRDS("SelectionDescriptionsA13A.rds")
# nPathways <- length(rownames(GenesER(rownames(seuratObj))@result))
#
# PrintOverlapCC(CCPathwaysOverlap(seuratObj, upDescriptions, selectionDescriptions,
#                                  RownamesCMLA13A, types, nPathways))
#
# View(upMarkers[[8]])
# View(selectionMarkers[[1]][intersect(rownames(upMarkers[[8]]), rownames(selectionMarkers[[1]])), ])
# View(upMarkers[[8]][intersect(rownames(selectionMarkers[[1]]), rownames(upMarkers[[8]])), ])
# View(selectionMarkers[[1]])
#
# #Section 10: Three-way overlaps
# x <- AllLitCCOverlap(geneSets$Lit, upMarkers, selectionMarkers,
#                      RownamesCMLA13A, types, length(rownames(seuratObj)))
# saveRDS(x, "ThreeWayLit.rds")
#
# x <- readRDS("ThreeWayLit.rds")
#
# x[[14]]
# PrintAllLitCCOverlap(x, LitNames)
# PrintOverlapCC(LogCCOverlap(geneSets$side, upMarkers, selectionMarkers, RownamesCMLA13A,
#                             types, length(rownames(seuratObj))))
#
# #Section 11: Epigenetics
# pathways <- readRDS("PathwaysA13A.rds")
# epiPathways <- FindEpigeneticPathways(pathways, allEpi)
# upEpiPathways <- FindClustersEpiPathways(upDescriptions, epiPathways)
# CountClusterEpiPathways(upEpiPathways, types)
# PrintOverlap(LitLogV(epiPathways, upDescriptions, types, length(pathways), F))
# clEpiPathways <- FindGroupingsEpiPathways(clusterDescriptions, epiPathways)
# CountGroupingEpiPathways(clEpiPathways, types)
# PrintGroupingEpiPathways(clEpiPathways, types)

# x <- LitMarkersOverlap(seuratObj, geneSets, upMarkers)
# PrintOverlap(x[[14]])
# SideMarkersOverlap(seuratObj, geneSets, upMarkers)
#
# cellMarker <- EnrichrVToER(enriched, "CellMarker_Augmented_2021")
# panglao <- EnrichrVToER(enriched, "PanglaoDB_Augmented_2021")
# ConceptCouple(cellMarker, 2, 6, 5, "CellMarker_Augmented_2021")
# ConceptCouple(panglao, 4, 7, 5, "PanglaoDB_Augmented_2021")
# jensen <- EnrichrVToER(enriched, "Jensen_DISEASES")
# StoreCnetplot(jensen[[11]])
#
# View(ConceptCouple)



