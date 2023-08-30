#' @export

#Obtain the desired ordering by using "unique" rather than "levels"
SplitWithRetainingZeros <- function(seuratObj)
  return(lapply(unique(seuratObj$orig.ident), function(x) subset(seuratObj, orig.ident == x)))

#Apply SCTransform individually to each Seurat subset corresponding to a treatment condition
IndividualSCT <- function(seuratObj){
  seurats <- SplitWithRetainingZeros(seuratObj)
  seurats <- lapply(seurats, function(seuratObj){
    seuratObj <- NormalizeData(seuratObj)
    seuratObj <- CellCycleScoring(seuratObj,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
    seuratObj <- SCTransform(seuratObj, vst.flavor="v2", vars.to.regress = c("Phase", "percent.mt", "percent.ribo"),
                             return.only.var.genes = F, min_cells = 0)
    return(seuratObj)
  })
  VarFeatList <- lapply(seurats, VariableFeatures)
  VarFeatUnion <- Reduce(union, VarFeatList)
  seurats <- PrepSCTIntegration(seurats, anchor.features = VarFeatUnion)
  seuratObj <- merge(seurats[[1]],seurats[c(2:length(seurats))])
  VariableFeatures(seuratObj) <- VarFeatUnion
  seuratObj <- PrepSCTFindMarkers(seuratObj)

  #Correcting NA introduced by PrepSCTFindMarkers
  seuratObj[["SCT"]]@counts[is.na(seuratObj[["SCT"]]@counts)] <- 0
  seuratObj[["SCT"]]@data[is.na(seuratObj[["SCT"]]@data)] <- 0

  #Removing rare genes induced by the SCT modification of the counts assay
  rareGenes <- FindRareGenes(seuratObj, assay = "SCT")
  seuratObj <- RemoveRareGenes(seuratObj, rareGenes)

  seuratObj$orig.ident <- as.factor(seuratObj$orig.ident)
  return(seuratObj)
}

IndividualNormscal <- function(seuratObj){
  seurats <- SplitObject(seuratObj, "orig.ident")
  seurats <- lapply(seurats, function(seuratObj){
    seuratObj <- NormalizeData(seuratObj)
    seuratObj <- FindVariableFeatures(seuratObj)
    seuratObj <- CellCycleScoring(seuratObj,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
    return(seuratObj)
  })
  VarFeatList <- lapply(seurats, VariableFeatures)
  VarFeatUnion <- Reduce(union, VarFeatList)
  seuratObj <- merge(seurats[[1]],seurats[c(2:length(seurats))])
  VariableFeatures(seuratObj) <- VarFeatUnion
  seuratObj <- ScaleData(seuratObj, vars.to.regress = c("orig.ident", "Phase", "nCount_RNA", "percent.mt"))

  seuratObj$orig.ident <- as.factor(seuratObj$orig.ident)
  return(seuratObj)
}
