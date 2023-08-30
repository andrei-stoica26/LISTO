#' @export

MarkersJaccard <- function(markers1, markers2, threshold, isLog = T)
  return(Jaccard(rownames(FilterMarkers(markers1, threshold, isLog)), rownames(FilterMarkers(markers2, threshold, isLog))))


TwoLogsJaccard <- function(markers1, markers2, N, isLog = T){
  if (length(rownames(markers1)) * length(rownames(markers2)) == 0)return(1)
  logs <- LogsMarkerLists(markers1, markers2, isLog = T)
  df <- data.frame(
    pvalue = unlist(lapply(logs, function(x) TwoLogsPV(markers1, markers2, x, N, isLog))),
    Jaccard = unlist(lapply(logs, function(x) MarkersJaccard(markers1, markers2, x, isLog))))
  df <- df[order(df$pvalue), ]
  df$pvalue <- BY(df$pvalue, 0.05)$Adjusted.pvalues
  return(df[median(1:length(rownames(df))), ])
}
