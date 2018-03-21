get_log_Pr_ck_matrix <- function(cl,X,Zck,r,S){
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

  # Add the BIC penalty
  BIC_penalty = (S_length * log(apply(Zck,2,sum)))/2
  log_Pr_ck_matrix = t(t(log_Pr_ck_matrix)+BIC_penalty)

  return(log_Pr_ck_matrix)
}
