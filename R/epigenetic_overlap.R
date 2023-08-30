#' @export

allEpi <- "histone|chromatin|epigenetic|protein-DNA|packaging|nucleosome|DNA methylation|imprinting|silencing|ncRNA-mediated|conformation|geometric|miRNA|inactivation of X chromosome"

FindPathways <- function(seuratObj, database = "enrichGO")
  return(GenesER(rownames(seuratObj), database)@result$Description)

FindEpigeneticPathways <- function(pathways, pattern)
  return (pathways[grep(pattern, pathways)])

FindGroupingsEpiPathways <- function(descriptions, epiPathways)
  return(lapply(descriptions, function(x) intersect(rownames(x), epiPathways)))

CountGroupingEpiPathways <- function(grEpiPathways, groupingNames){
  df <- data.frame(Count = sapply(1:length(groupingNames), function(i) length(grEpiPathways[[i]])))
  rownames(df) <- groupingNames
  df <- subset(df, Count > 0)
  df <- df[order(df$Count, decreasing = T),, drop = F]
  return(df)
}

PrintGroupingEpiPathways <- function(grEpiPathways, groupingNames)
  return(sapply(1:length(groupingNames), function(i){
    print(groupingNames[i])
    print(grEpiPathways[[i]])
  }))

EpigeneticSolo <- function(enrichments, index, epiPathways, types, typeClass, cutoff = NA){
  aux <- enrichments[[index]]
  aux@result <- subset(aux@result, Description %in% epiPathways)
  dev.new(width =  8, height = 8, noRStudioGD = TRUE)
  if (is.na(cutoff))
    return(StoreCnetplot(aux, length(aux@result$Description)) + ggtitle(str_c("Epigenetic GO terms enriched in the ",
                                                                              types[index], " ", typeClass)) +
             theme(plot.title = element_text(hjust = 0.5, size = 10)))
  return(StoreCnetplot(aux, cutoff) + ggtitle(str_c("Top ", cutoff, " epigenetic GO terms enriched in the ",
                                                    types[index], " ", typeClass)) +
           theme(plot.title = element_text(hjust = 0.5, size = 10)))
}

EpigeneticCnetplot <- function(enrichments, index1, index2, epiPathways, title){
  aux1 <- enrichments[[index1]]
  aux1@result <- subset(enrichments[[index1]]@result, Description %in% epiPathways)
  aux2 <- enrichments[[index2]]
  aux2@result <- subset(enrichments[[index2]]@result, Description %in% epiPathways)
  return(GrobPlot(list(StoreCnetplot(aux1, length(aux1@result$Description)) + NoLegend(),
                       StoreCnetplot(aux2, length(aux2@result$Description))+ NoLegend()), 2, leftMargin = 5) + ggtitle(title) +
           theme(plot.title = element_text(hjust = 0.5)))
}
