split_clusters <- function(X,Zck,r){
  # For each cluster
  for(clusterID in dim(Zck)[2]){
    # Get all cells in the cluster
    Xc=X[,Zck[,clusterID]]
    # First, calculate deltaYgt for all g (genes) and t (thresholds)
    library(plyr)
    tmp2 = ldply(parLapply(cl,as.list(1:dim(Xc)[1]),FUN=get.deltaYgt,Xc,r),data.frame)
    # - Need to find optimal values of t for each gene
    # - Find the ten genes giving top values of deltaYgt
    # For each of the top ten genes, split
  }
}
