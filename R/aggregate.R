#' @export

#Functions for the aggregate analysis of two Seurat datasets

SeuratCoupleSynonyms <- function(seuratObj1, seuratObj2){
  questGenes <- setdiff(rownames(seuratObj1), rownames(seuratObj2))
  synonyms1 <- c()
  synonyms2 <- c()
  for (x in questGenes){
    #Avoiding duplicates, hence the setdiff
    foundGene <- setdiff(intersect(unlist(humanSyno(x)), rownames(seuratObj2)), rownames(seuratObj1))
    if (length(foundGene) == 1){
      synonyms1 <- c(synonyms1, foundGene)
      synonyms2 <- c(synonyms2, x)
    }
  }
  df <- data.frame(foundGene = synonyms2)
  rownames(df) <- synonyms1
  return(df)
}

SynonymizeNames <- function(geneNames, synonymsDF){
  for (i in 1:length(geneNames)){
    if (geneNames[i] %in% rownames(synonymsDF))
      geneNames[i] = synonymsDF[geneNames[i], ]
  }
  return(geneNames)
}

SynonymizeMarkers <- function(markerList, synonymsDF)
  return(lapply(markerList, function(x){
    rownames(x) <- SynonymizeNames(rownames(x), synonymsDF)
    return(x)
  }))
