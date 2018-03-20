calculate_cluster_mean <- function(X,Zck,cl){
  # Error check
  if(class(Zck[,1])!="logical"){stop("Zck must be logical, not numeric")}

  get.mean.geneexp <- function(whichCells,X){
    regularised_mean <- function(x){
      out = (10^-4 + sum(x)) / (1+length(x))
    }
    #print(which(whichCells))
    #print(X[,whichCells,drop=FALSE])
    return(apply(X[,whichCells,drop=FALSE],1,regularised_mean))
  }
  Ugk = parApply(cl,Zck,2,get.mean.geneexp,X)

  return(Ugk)
}
