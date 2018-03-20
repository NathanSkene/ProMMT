# Expectation step
expectation_step <- function(cl,X,Zck,r,S){
  # z_c_k = pr(k|X_c,S,u)
  Ugk  = calculate_cluster_mean(X,Zck,cl)
  pi_k = calculate_priors(Zck) # Set the priors
  Pgk  = calculate_Pgk(Ugk,r)

  # Calculate log Probability of X_c belonging to cluster K
  get_logPr <- function(params,X,Pgk,r,pi_k,S){
    k=params["k"]
    c=params["c"]
    logSumFunc <- function(g,X,Pgk,r){
      return(X[g,c]*log(Pgk[g,k]) + r*log(1-Pgk[g,k]))
    }
    return(log(pi_k[k]) + sum(unlist(lapply(as.list(S),FUN=logSumFunc,X,Pgk,r))))
  }

  paramList = list()
  for(k in 1:dim(Zck)[2]){for(c in 1:dim(X)[2]){paramList[[length(paramList)+1]]=c(k=k,c=c)}}
  logPr_ck = parLapply(cl=cl,paramList,get_logPr,X,Pgk,r,pi_k,S)
  log_Pr_ck_matrix = matrix(unlist(logPr_ck),nrow=dim(X)[2],ncol=dim(Zck)[2])

  Zck_new = get_best_clust_given_logPr_ck(log_Pr_ck_matrix)

  # If a cluster now has no cells, then drop that cluster
  Zck_new = Zck_new[,!colSums(Zck_new)==0]

  # Force Zck to be logical
  Zck_new = Zck_new>0

  return(list(Zck=Zck_new,log_Pr_ck_matrix=log_Pr_ck_matrix))
}
