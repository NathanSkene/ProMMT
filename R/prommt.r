# Run the ProMTT algorithm
ProMTT <- function(X,cl){
  # Initial setup
  S_length = 150 # How many genes shoul be permitted to have different means?
  S = sample(1:dim(X)[1],S_length) # Numerical index of genes whose means are permitted to vary
  r=2 # Value of r (for negative binomial parameter) is fixed
  Zck  = diag(dim(X)[2])>0 # Each cell is assigned to it's own class

  noChange=0
  clusters = list()
  for(iter in 0:200){
    #if(iter!=0){
    #  Zck = split_clusters
    #}

    # EM steps
    expOut = ProMMT::expectation_step(cl,X,Zck,r,S)
    dim(expOut$Zck)
    print(colSums(expOut$Zck))
    Zck = expOut$Zck
    S_old = S
    S   = ProMMT::maximisation_step(cl,X,Zck,r,S,S_length)
    print(sprintf("Percent of S replaced: %s",100-(length(intersect(S,S_old))/length(S_old))*100))

    # Calculate BIC penalty
    BIC_penalty = (S_length * log(apply(Zck,2,sum)))/2

    # Prune
    prune_out = ProMMT::prune_clusters(X,Zck,r,S,BIC_penalty)

    if(prune_out$dropped!=-1){
      Zck_old = Zck
      Zck = prune_out$Zck
      print(sprintf("------OUTPUT: penalty: %s, diff: %s",prune_out$BIC_penalty,prune_out$diffLogLik))
    }

    # Get score for this model & add to list
    log_Pr_ck_matrix = ProMMT::get_log_Pr_ck_matrix(cl,X,Zck,r,S)
    log_Pr = sum(log_Pr_ck_matrix[Zck])
    clusters[[length(clusters)+1]] = list(Zck=Zck,log_Pr=log_Pr,log_Pr_ck_matrix=log_Pr_ck_matrix)

    # Should loop be terminated? (have clusters stopped changing?)
    if(prune_out$dropped==-1){
      noChange=noChange+1
      if(noChange>5){break()}
    }else{
      noChange=0
    }
  }

  logPrs = unlist(lapply(clusters,function(x) x$log_Pr))
  plot(1:length(logPrs),logPrs)
  useClust = which(logPrs==max(logPrs))[1]
  out = list()
  out$top_clust = clusters[[useClust]]
  out$clusters  = clusters
  out$logPrs    = logPrs
  return(out)
}
