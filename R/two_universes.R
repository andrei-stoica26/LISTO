#' @export

#Functions for overlap assessments for two scRNA-seq datasets
Choose <- function(n, k)
  return(Legendre(n) - (Legendre(k) + Legendre(n-k)))

Legendre <- function(n){
  result <- rep(0, N_PRIMES)
  for (i in 1:N_PRIMES){
    k <- PRIMES[i]
    while (k <= n){
      result[i] <- result[i] + floor(n/k)
      k <- k * PRIMES[i]
    }
  }
  return(result)
}


PowerProduct <- function(exponents)
  return(prod(mapply(function(x, y)x^y, PRIMES, exponents)))

MNDenom <- function(m_int_n, b_int_m)
  return(Choose(m_int_n, b_int_m))


MNCount <- function(a_int_n, a_int_b, m_int_n, b_int_m, denom = MNDenom(m_int_n, b_int_m)){
  exponents <- Choose(a_int_n, a_int_b) + Choose(m_int_n - a_int_n, b_int_m - a_int_b) - denom
  return(PowerProduct(exponents))
}

MNPVal <- function(a_int_n, a_int_b, m_int_n, b_int_m){
  denom <- MNDenom(m_int_n, b_int_m)
  return(sum(unlist(lapply(a_int_b:min(a_int_n, b_int_m), function(i)MNCount(a_int_n, i, m_int_n, b_int_m, denom)))))
}

MNGenesCount <- function(M, N, a, b)
  return(MNCount(length(intersect(a, N)),
                 length(intersect(a, b)),
                 length(intersect(M, N)),
                 length(intersect(b, M))
  ))

MNGenesPVal <- function(M, N, a, b)
  return(MNPVal(length(intersect(a, N)),
                length(intersect(a, b)),
                length(intersect(M, N)),
                length(intersect(b, M))
  ))

MNTwoLogsPV <- function(M, N, markers1, markers2, threshold, isLog = T)
  return(MNGenesPVal(M, N, rownames(FilterMarkers(markers1, threshold, isLog)), rownames(FilterMarkers(markers2, threshold, isLog))))

MNTwoLogs <- function(M, N, markers1, markers2, isLog = T){
  if (length(rownames(markers1)) * length(rownames(markers2)) == 0)return(1)
  logs <- LogsMarkerLists(markers1, markers2, isLog)
  return(median(BY(unlist(lapply(logs, function(x) MNTwoLogsPV(M, N, markers1, markers2, x, isLog))), 0.05)$Adjusted.pvalues))
}
MNTwoLogsParallel <- function(clust, M, N, markers1, markers2, isLog = T){
  if (length(rownames(markers1)) * length(rownames(markers2)) == 0)return(1)
  logs <- LogsMarkerLists(markers1, markers2, isLog)
  return(median(BY(unlist(parLapply(clust, logs, function(x) MNTwoLogsPV(M, N, markers1, markers2, x, isLog))), 0.05)$Adjusted.pvalues))
}


MNOverlap <- function(markers1, markers2, typeNames1, typeNames2, M, N, isLog = T){
  pvalues <- c()
  pairs <- c()
  for (i in 1:length(markers1))
  {
    message(str_c("Calculating overlaps with cluster ", i - 1, " from A13A Seurat object"))
    for (j in 1:length(markers2)){
      message(str_c("Calculating overlaps with cluster ", j - 1, " from patient Seurat object"))
      pairs <- c(pairs, str_c(typeNames1[i], " (A13A) and ", typeNames2[j],  " (Patient)"))
      pvalues <- c(pvalues, MNTwoLogsParallel(clust, M, N, markers1[[i]], markers2[[j]], isLog))
    }
  }
  df <- data.frame(Grouping = pairs, pvalue = pvalues)
  return(BYCorrectDF(df))
}

FixTU <- function(markersTU){
  df <- as.data.frame(do.call(rbind, str_split(markersTU$Grouping, " and ")))
  df$pvalue <- markersTU$pvalue
  colnames(df) <- c("Cluster", "Grouping", "pvalue")
  return(df)
}
