#' @export

AddDiversity <- function(seuratObj, type, typeStr){
  countsMat = as.matrix(GetAssayData(seuratObj, "counts"))
  currentDiversity <- vegan::diversity(countsMat, type, MARGIN = 2)
  message(str_c("Calculated ", typeStr, " diversity."))
  seuratObj@meta.data[[str_c(type, ".diversity")]] <- currentDiversity
  message (str_c("Added ", typeStr, " diversity to Seurat object."))
  gc()
  return (seuratObj)
}

AddShannonDiversity <- function(seuratObj)
  AddDiversity(seuratObj, "shannon", "Shannon")

AddSimpsonDiversity <- function(seuratObj)
  AddDiversity(seuratObj, "simpson", "Simpson")

AddFilteringCriteria <- function(seuratObj, manyCriteria = T){
  if (manyCriteria){
    seuratObj <- AddShannonDiversity(seuratObj)
    seuratObj <- AddSimpsonDiversity(seuratObj)
    seuratObj <- PercentageFeatureSet(seuratObj, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")
    seuratObj$novelty <- log10(seuratObj$nFeature_RNA) / log10(seuratObj$nCount_RNA)
  }
  seuratObj <- PercentageFeatureSet(seuratObj, pattern = "^MT-", col.name = "percent.mt")
  return (seuratObj)
}

#Preventing errors due to empty Seurat subsets
SafeSeurat <- function(seuratObj)
  tryCatch(seuratObj, error = function(...) NULL)

PruneSeurat <- function(seuratObj, criteria, msg){
  removedSeurat <- SafeSeurat(seuratObj[, which(criteria)])
  message(msg)
  CountCells(removedSeurat, "out")
  keptSeurat <- SafeSeurat(seuratObj[, which(!criteria)])
  CountCells(keptSeurat, "in")
  return(keptSeurat)
}

FetchFilter <- function (seuratObj, column, lowerBound, upperBound, text){
  column <- FetchData(seuratObj, vars = column)
  seuratObj <- PruneSeurat(seuratObj, column < lowerBound | column > upperBound,
                           str_c("Filtering the Seurat object: ", text))
  return(seuratObj)
}

GateCells <- function(seurats, colStr, bounds, texts)
  if (!is.null(bounds))
    mapply(function(x, y, z, t)
      FetchFilter(x, colStr, y[[1]], y[[2]], str_c(z, " per cell - ", t, ".")), seurats, bounds, texts, conditions)

QualityControl <- function(seuratObj, doRemoveDoublets = F, noveltyBounds = NULL, mtBounds = NULL, riboBounds = NULL,
                           countBounds = NULL, featureBounds = NULL, shannonBounds = NULL, simpsonBounds = NULL){
  if(doRemoveDoublets)seuratObj <- PruneSeurat(seuratObj, FetchData(seuratObj, vars = "unit.class") == "doublet",
                                               "Removing doublets.")
  seurats <- SplitObject(seuratObj, "orig.ident")
  allCols <- c("novelty", "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "shannon.diversity", "simpson.diversity")
  allBounds <- list(noveltyBounds, featureBounds, countBounds, mtBounds, riboBounds, shannonBounds, simpsonBounds)
  allTexts <- c("Novelty", "Number of genes", "Number of UMIs", "Percentage of mitochondrial genes", "Percentage of ribosomal genes",
                "Shannon diversity", "Simpson diversity")
  for (i in 1:7) seurats <- GateCells(seurats, allCols[i], allBounds[[i]], allTexts[i])
  seuratObj <- merge(seurats[[1]],seurats[c(2:length(seurats))])
  message("Finished.")
  return(seuratObj)
}


FindRareGenes <- function(seuratObj, assay = "RNA", nCells = 10){
  expMatrix <- seuratObj[[assay]]@counts
  df <- data.frame(Count = rowSums(expMatrix != 0))
  df <- subset(df, Count < nCells)
  return(df)
}

RemoveRareGenes <- function(seuratObj, rareGenes){
  retainedGenes <- setdiff(rownames(seuratObj), rownames(rareGenes))
  seuratObj <- subset(seuratObj, features = retainedGenes)
  message(str_c(length(rownames(rareGenes)), " rare genes removed."))
  message(str_c(length(retainedGenes), " genes retained in the Seurat object."))
  return(seuratObj)
}
