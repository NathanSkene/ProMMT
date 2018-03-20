# Maximisation step
maximisation_step <- function(cl,X,Zck,r,S,S_length){
  # S and u are optimised given current values of z_c_k
  Ugk  = calculate_cluster_mean(X,Zck,cl)
  Ug0  = apply(X,1,mean)
  Pg0  = Ug0/(r+Ug0)
  Pgk  = calculate_Pgk(Ugk,r)

  # Calculate for each gene, Yg (the gain in log likelihood when the distribution of gene g is allowed to vary between classes)
  calc_Yg <- function(g,Zck,X,Pgk,Pg0,r){
    val = 0
    for(c in 1:dim(Zck)[2]){
      val = val + X[g,c]*(log(Pgk[g,c] - log(Pg0[g]))) + r*(log(1-Pgk[g,c]) - log(1-Pg0[g]))
    }
    return(val)
  }
  Yg = unlist(parLapply(cl,as.list(1:dim(X)[1]),FUN=calc_Yg,Zck,X,Pgk,Pg0,r))

  # Select the new set of S (those with biggest values of Yg)
  top_Yg_genes = names(sort(Yg,decreasing=TRUE)[1:S_length])
  S_new = which(rownames(X) %in% top_Yg_genes)
  return(S_new)
}
