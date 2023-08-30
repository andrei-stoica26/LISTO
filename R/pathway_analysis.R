#' @export

#Different databases are accessed in slightly different ways
GetEnrichmentResult <- function(entrezList, funString)
  switch(funString,
         "enrichDGN" = setReadable(enrichDGN(entrezList), 'org.Hs.eg.db', 'ENTREZID'),
         "enrichGO" = enrichGO(entrezList, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE),
         "enrichWP" = setReadable(enrichWP(entrezList, "Homo sapiens"), 'org.Hs.eg.db', 'ENTREZID'),
         "enrichKEGG" = setReadable(enrichKEGG(entrezList, "hsa"), 'org.Hs.eg.db', 'ENTREZID'),
         "enrichNCG" = setReadable(enrichNCG(entrezList), 'org.Hs.eg.db', 'ENTREZID'),
  )

EntrezGenes <- function(markerNames)
  AnnotationDbi::select(org.Hs.eg.db, keys = markerNames,
                        columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")[[2]]

EnrichedTerm <- function(funString){
  return (switch(funString,
                 "enrichDGN" = "Enriched disease term",
                 "enrichGO" = "Enriched GO term",
                 "enrichWP" = "Enriched WikiPathways term",
                 "enrichKEGG" = "Enriched KEGG term"
  )
  )
}

GenesER <- function(genes, database = "enrichGO")
  GetEnrichmentResult(EntrezGenes(genes), database)

StoreCnetplot <- function(enrichmentResult, nCategories = 20, database = "enrichGO", nodeLabel = "all")
  cnetplot(enrichmentResult, showCategory = nCategories, colorEdge = T, cex_label_category = 0.5,
           cex_label_gene = 0.5, node_label = nodeLabel) +
  scale_colour_manual("type", values = c("red", "purple")) + NoLegend()

#Plotting two concept network plots combined
ConceptCouple <- function(clusterEnrichments, i, j, nPathways, database){
  dev.new(width =  12, height = 12, noRStudioGD = TRUE)
  return(GrobPlot(list(StoreCnetplot(clusterEnrichments[[i + 1]], nPathways),
                       StoreCnetplot(clusterEnrichments[[j + 1]], nPathways)), 2, leftMargin = 5) +
           ggtitle(str_c(database, " concept network plot of markers of Cluster ", i, " (A) and Cluster ", j, " (B)")) +
           theme(plot.title = element_text(hjust = 0.5))
  )
}

#Plotting one concept network plot
ConceptSolo <- function(clusterEnrichments, i, nPathways, database){
  dev.new(width =  6, height = 6, noRStudioGD = TRUE)
  return(StoreCnetplot(clusterEnrichments[[i + 1]], nPathways) + ggtitle(str_c(database, " concept network plot of markers of Cluster ", i)) +
           theme(plot.title = element_text(hjust = 0.5, size = 10)))
}
