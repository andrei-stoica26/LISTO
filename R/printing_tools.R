#' @export

#Easy printing of gene lists and overlap matrices

MessageVector <- function(vect)
  invisible(lapply(vect, message))

PrintVector <- function(vect)
  invisible(lapply(vect, print))

PrintOverlap <- function(df){
  print(df)
  print(df$pvalue)
  MessageVector(df$Grouping)
}

PrintOverlap2 <- function(df){
  print(df)
  PrintVector(df$pvalue)
  MessageVector(df$Grouping)
}

PrintOverlapCC <- function(df){
  print(df)
  MessageVector(df$Cluster)
  PrintOverlap2(df)
}
