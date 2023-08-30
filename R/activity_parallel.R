#' @export

#Parallelizing the calculation of ORIGINS activity scores

ColumnSlicesDF <- function(df, nSlices){
  #Divides a dataframe by columns into roughly equally sized slices
  sliceSize <- floor(ncol(df)/nSlices)
  slices <- lapply(1:(nSlices - 1), function(x)df[, ((x - 1)*sliceSize + 1):(x * sliceSize)])
  slices <- append(slices, list(df[, ((nSlices - 1) * sliceSize + 1):ncol(df)]))
  return(slices)
}

RawActivitySingleCell <- function(singleCell, adjacency_list)
  return(sum(unlist(lapply(1:nrow(adjacency_list), function(gene){
    gene_name <- adjacency_list[gene, 1]
    return(singleCell[gene_name] * sum(singleCell[unlist(adjacency_list[gene_name, 2])]))}))))

ActivitySliceRun <- function(clust, dfSlice, nSlice)
  return(unlist(parLapply(clust, 1:ncol(dfSlice), function(x){
    singleCell <- dfSlice[,x]
    names(singleCell) <- rownames(dfSlice)
    RawActivitySingleCell(singleCell, adjacency_list)
  })))

MinMaxNormalization <- function(v)
  return((v - min(v)) / (max(v) - min(v)))

ActivityParallel <- function(expression_matrix, adjacency_edges = .data$differentiation_edges, nSlices = 20, nCores = 8){
  #Preparation for the ORIGINS activity run
  adjacency_edges <- differentiation_edges
  cc.df <- data.frame(adjacency_edges)
  colnames(cc.df) <- c("node1", "node2")
  cc.df <- cc.df %>% dplyr::arrange((.data$node1))
  reduced.df <- data.frame(expression_matrix)
  reduced.df <- reduced.df %>% dplyr::arrange((row.names(reduced.df)))
  intersection_genes <- intersect(rownames(reduced.df), cc.df$node1)
  cc.df_reduced <- cc.df[cc.df$node1 %in% intersection_genes, ]
  cc.df_reduced <- cc.df_reduced[cc.df_reduced$node2 %in% intersection_genes, ]
  adjacency_list <- stats::aggregate(node1 ~ node2, cc.df_reduced, FUN = I)
  rownames(adjacency_list) <- adjacency_list[, 1]
  dfSlices <- ColumnSlicesDF(reduced.df, nSlices)
  clust <- makeCluster(nCores)
  clusterExport(clust, c("RawActivitySingleCell", "adjacency_list"),
                envir=environment())
  activitySlices <- lapply(1:nSlices, function(x){
    dfSlice <- dfSlices[[x]]
    clusterExport(clust, "dfSlice",envir=environment())
    start_time <- Sys.time()
    activitySlice <- ActivitySliceRun(clust, dfSlice)
    end_time <- Sys.time()
    timeElapsed <- end_time - start_time
    message(str_c("Slice ", x, " finished in ", timeElapsed, " seconds."))
    return(activitySlice)
  })
  stopCluster(clust)
  act <- unlist(activitySlices)
  return(MinMaxNormalization(act))
}
