#' @export

#Fitting a gam using the gam function

FitGam <- function(seuratObj, expression, genes, column)
  return(apply(expression[genes, ], 1, function(z)
    return(summary(gam(z ~ lo(seuratObj[[column]]),
                       data = data.frame(z = z, t = seuratObj[[column]])))[4][[1]][1,5])))

SaveGam <- function(gam, fileName)
  saveRDS(data.frame(pvalue = sort(gam)), fileName)

PrepareGamForOverlapTests <- function(gam){
  #Adjust p-values with Bonferroni and subset
  gam$pvalue <- gam$pvalue * length(rownames(gam))
  gam <- subset(gam, pvalue < 0.05)
  #label just makes it compatible with the column name normally used in overlap assessments
  gam$avg_log2FC <- -log10(gam$pvalue)
  return(gam)
}

GetGeneIndices <- function(seuratObj, df)
  return(sapply(rownames(seuratObj), function(x) which(rownames(df) %in%  x)))



