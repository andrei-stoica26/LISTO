#' @export

#Generating overlap matrices
#A V in a function name means the corresponding overlap, denoted by the preceding elemenet, involved a list
#For example: LitLogV computes the overlap between a single set of genes and multiple sets of markers (e.g. cluster markers)
#But LitVLogV computes the overlap between multiple sets of genes and multiple sets of markers (e.g. cluster markers)
#BYCorrect performs a Benjamini-Yekutieli correction

LitLogV <- function(litMarkers, logVMarkers, groupingNames, N, isLog = T){
  df <- data.frame(Grouping = groupingNames, pvalue = unlist(lapply(logVMarkers,
                                                                    function(x)OneLitOneLog(litMarkers, x, N, isLog))))
  return (BYCorrectDF(df))
}

LitVLogV <- function(litVMarkers, logVMarkers, groupingNames, N, isLog = T)
  return(lapply(litVMarkers, function(x) LitLogV(x, logVMarkers, groupingNames, N, isLog)))

LogLogV <- function(logMarkers, logVMarkers, groupingNames, N, isLog = T){
  df <- data.frame(Grouping = groupingNames, pvalue = unlist(lapply(logVMarkers,
                                                                    function(x)TwoLogs(logMarkers, x, N))))
  return (BYCorrectDF(df))
}

LitMarkersOverlap <- function(seuratObj, geneSets, markerList, groupingNames)
  return(LitVLogV(geneSets$Lit, markerList, groupingNames, length(rownames(seuratObj))))

SideMarkersOverlap <- function(seuratObj, geneSets, markerList, groupingNames)
  return(LogLogV(geneSets$side, markerList, groupingNames, length(rownames(seuratObj))))

LitLitV <- function(litMarkers, litVMarkers, groupingNames, N){
  df <- data.frame(Grouping = groupingNames, pvalue = unlist(lapply(litVMarkers,
                                                                    function(x)TwoGeneSetsPVal(litMarkers, x, N))))
  return (BYCorrectDF(df))
}

