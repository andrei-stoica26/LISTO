#' @export

#Convert Enrichr output to enrichment result type of object, for plotting purposes

EnrichrToER <- function(enrichr, nTests){
  enrichr$Genes <- sapply(enrichr$Genes, function(x) gsub(";", "/", x))
  #Correcting with Bonferroni for multiple testing
  enrichr$Adjusted.P.value <- enrichr$Adjusted.P.value * nTests
  enrichr <- subset(enrichr, Adjusted.P.value < 0.05)
  if (length(rownames(enrichr)) > 0)
    return(enrichDF2enrichResult(enrichr,
                                 keyColname = "Term",
                                 descriptionColname = "Term",
                                 pvalueColName = "Adjusted.P.value",
                                 geneColName = "Genes",
    )) else return (enrichr)
}

EnrichrVToER <- function(enrichrV, database)
  return(lapply(enrichrV, function(x) EnrichrToER(x[[database]], length(enrichrV))))

