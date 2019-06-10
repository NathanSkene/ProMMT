split_clusters <- function(cl,X,S_length,Zck,r){
  # For each cluster
  for(clusterID in dim(Zck)[2]){
    # Get all cells in the cluster
    Xc=X[,Zck[,clusterID]]
    # First, calculate deltaYgt for all g (genes) and t (thresholds)
    library(plyr)
    library(magrittr)
    deltaYgt = ldply(parLapply(cl,as.list(1:dim(Xc)[1]),FUN=get.deltaYgt,Xc,r),data.frame) %>%
      arrange(desc(deltaYgt))
    # - Need to find optimal values of t for each gene
    tenGenesWithBiggestChange = deltaYgt %>% .[!duplicated(.$g),] %>% .[1:10,]
    # - Find the ten genes giving top values of deltaYgt
    # For each of the top ten genes, split
    inputGeneThresh = tenGenesWithBiggestChange[1,]
    splitZck = data.frame(clust1=Xc[inputGeneThresh$g,]>=inputGeneThresh$t,clust2=Xc[inputGeneThresh$g,]<inputGeneThresh$t)
    converged = ProMMT::ProMTT(cl=cl,X=X,S_length=S_length,r=r,Zck=splitZck)
  }
}
