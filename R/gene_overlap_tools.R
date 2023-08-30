#' @export

#Functions for calculating overlaps between gene sets.
#Lists that serve as input can use logs for ranking, in which case isLog must be set to TRUE.
#They can also use p-values, in which case isLog must be set to FALSE.

################### Utilities
FilterMarkers <- function(markers, cutoff, isLog = T)
  if (isLog) return (subset(markers, markers$avg_log2FC > cutoff)) else
    return (subset(markers, markers$avg_log2FC < cutoff))


################### Two gene sets
TwoSetsPVal <- function(x, y, N, lowerTail){
  lx <- length(x)
  ly <- length(y)
  lint <- length(intersect(x, y))
  return (phyper(lint - (1 - lowerTail), lx, N - lx, ly, lower.tail = lowerTail))
}

#Gene sets: we look at the probability that the gene sets intersect in at least intersect(x, y) points
TwoGeneSetsPVal <- function(x, y, N)
  return(TwoSetsPVal(x, y, N, F))



TwoCellSetsPVal <- function(x, y, N)
  return(TwoSetsPVal(x, y, N, T))
#Gene sets: we look at the probability that the gene sets intersect in at most intersect(x, y) points


OneLitOneLogPV <- function(litMarkers, markers, threshold, N, isLog = T)
  return(TwoGeneSetsPVal(litMarkers, rownames(FilterMarkers(markers, threshold, isLog)), N))

OneLitOneLog <- function(litMarkers, markers, N, isLog = T){
  if (length(rownames(markers)) == 0)return(1)
  logs <- unique(sort(c(1 - isLog, markers$avg_log2FC)))
  return(median(BY(unlist(lapply(logs, function(x) OneLitOneLogPV(litMarkers, markers, x, N, isLog))), 0.05)$Adjusted.pvalues))
}

LogsMarkerLists <- function(markers1, markers2, isLog = T){
  logs <- unique(sort(c(1 - isLog, markers1$avg_log2FC, markers2$avg_log2FC)))
  if (isLog)
    logs <- subset(logs, logs < min(max(markers1$avg_log2FC), max(markers2$avg_log2FC))) else
      logs <- subset(logs, logs > max(min(markers1$avg_log2FC), min(markers2$avg_log2FC)))
    return (logs)
}

TwoLogsPV <- function(markers1, markers2, threshold, N, isLog = T)
  return(TwoGeneSetsPVal(rownames(FilterMarkers(markers1, threshold, isLog)), rownames(FilterMarkers(markers2, threshold, isLog)), N))

TwoLogs <- function(markers1, markers2, N, isLog = T){
  if (length(rownames(markers1)) * length(rownames(markers2)) == 0)return(1)
  logs <- LogsMarkerLists(markers1, markers2, isLog = T)
  return(median(BY(unlist(lapply(logs, function(x) TwoLogsPV(markers1, markers2, x, N, isLog))), 0.05)$Adjusted.pvalues))
}

ThreeSetsExact <- function(a, b, c, k, N){
  #The probability that three sets of cardinalities a, b and c intersect in exactly k points
  exact <- sum(sapply(max(a + b - N, k):min(a, b, N + k - c),
                      function(x) dhyper(x, a, N - a, b) * dhyper(k, x, N - x, c)))
  return (exact)
}

TestThreeSetsExact <- function(a, b, c, N)
  #Must return 1
  return(sum(sapply(0:min(a,b,c), function(x) ThreeSetsExact(a, b, c, x, N))))

ThreeGeneSetsExact <- function(smallSet, midSet, largeSet, k, N){
  lengths <- sort(c(length(smallSet), length(midSet), length(largeSet)))
  return (ThreeSetsExact(lengths[[1]], lengths[[2]], lengths[[3]], k, N))
}

ThreeGeneSetsPVal <- function(smallSet, midSet, largeSet, N){
  lowerTailBound <- length(Reduce(intersect, list(smallSet, midSet, largeSet)))
  if (lowerTailBound == 0) return(1)
  lengths <- sort(c(length(smallSet), length(midSet), length(largeSet)))
  #It takes the upper tail directly, which is more computationally expensive than taking the lower tail and subtracting it from 1
  #But the latter approach leads to precision errors
  return (sum(sapply(lowerTailBound:lengths[[1]], function(x) ThreeSetsExact(lengths[[1]], lengths[[2]], lengths[[3]], x, N))))
}

OneLitTwoLogsPV <- function(litMarkers, markers1, markers2, threshold, N, isLog = T)
  return(ThreeGeneSetsPVal(litMarkers, rownames(FilterMarkers(markers1, threshold)), rownames(FilterMarkers(markers2, threshold, isLog)), N))

OneLitTwoLogs <- function(litMarkers, markers1, markers2, N, isLog = T){
  logs <- LogsMarkerLists(markers1, markers2)
  return(median(BY(unlist(lapply(logs, function(x) OneLitTwoLogsPV(litMarkers, markers1, markers2, x, N, isLog))), 0.05)$Adjusted.pvalues))
}

LogsThreeMarkerLists <- function(markers1, markers2, markers3, isLog = T){
  logs <- sort(unique(c(1 - isLog, markers1$avg_log2FC, markers2$avg_log2FC, markers3$avg_log2FC)))
  if (isLog)
    logs <- subset(logs, logs < min(max(markers1$avg_log2FC), max(markers2$avg_log2FC), max(markers3$avg_log2FC))) else
      logs <- subset(logs, logs > max(min(markers1$avg_log2FC), min(markers2$avg_log2FC), min(markers3$avg_log2FC)))
    return (logs)
}

ThreeLogsPV <- function(markers1, markers2, markers3, threshold, N, isLog = T)
  return(ThreeGeneSetsPVal(rownames(FilterMarkers(markers1, threshold)), rownames(FilterMarkers(markers2, threshold, isLog)),
                           rownames(FilterMarkers(markers3, threshold)), N))

ThreeLogs <- function(markers1, markers2, markers3, N, isLog = T){
  logs <- LogsThreeMarkerLists(markers1, markers2, markers3)
  return(median(BY(unlist(lapply(logs, function(x) ThreeLogsPV(markers1, markers2, markers3, x, N, isLog = T))), 0.05)$Adjusted.pvalues))
}

BYCorrectDF <- function(overlapDF){
  overlapDF <- overlapDF[order(overlapDF$pvalue), ]
  overlapDF$pvalue <- BY(overlapDF$pvalue, 0.05)$Adjusted.pvalues
  overlapDF <- subset(overlapDF, pvalue < 0.05)
  return(overlapDF)
}
