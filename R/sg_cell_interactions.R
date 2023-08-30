#' @export

#SingleCellSignalR

MVInteractions <- function(data,genes,cluster,c.names=NULL,n=30){
  if (is.null(c.names)==TRUE){
    c.names <- paste("cluster",seq_len(max(cluster)))
  }
  rownames(data) <- genes
  l_sc <- LRdb[is.element(LRdb$ligand,rownames(data)),]
  int_sc <- l_sc[is.element(l_sc$receptor,rownames(data)),]
  lr_sc <- matrix(0,nrow=nrow(int_sc),ncol=max(cluster)^2)
  rownames(lr_sc) <- paste(int_sc$ligand,int_sc$receptor,sep=" / ")
  med <- sum(data)/(nrow(data)*ncol(data))
  nam <- NULL
  q <- 0
  for (i in seq_len(max(cluster))){
    for (j in seq_len(max(cluster))){
      q <- q+1
      lr_sc[,q] <- (rowMeans(data[int_sc$ligand,cluster==i])*
                      rowMeans(data[int_sc$receptor,cluster==j]))^0.5/
        (med + (rowMeans(data[int_sc$ligand,cluster==i])*rowMeans(
          data[int_sc$receptor,cluster==j]))^0.5)
      nam=c(nam,paste(c.names[i],c.names[j],sep=" -> "))
    }
  }
  colnames(lr_sc) <- nam
  if (sum(lr_sc)!=0){
    if (nrow(lr_sc)<n){
      n <- nrow(lr_sc)
    }
    lr_sc <- subset(lr_sc,rowSums(lr_sc)!=0)
    v <- apply(lr_sc,1,var)/apply(lr_sc,1,mean)
    lr_sc <- lr_sc[order(v,decreasing=TRUE),]
    lr_sc <- lr_sc[apply(lr_sc,1, max)>0.5,]
    pheatmap(t(lr_sc[seq_len(n),colSums(lr_sc[seq_len(n),])!=0]),
             cluster_cols=TRUE, fontsize = 7, treeheight_col = 0)
  } else {
    cat("No interactions detected. Make sure the genes vector is composed of
        HUGO official gene names.",fill=TRUE)
  }
}
