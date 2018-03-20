# Run the ProMTT algorithm
ProMTT <- function(X,cl){
  # Initial setup
  S_length = 150 # How many genes shoul be permitted to have different means?
  S = sample(1:dim(X)[1],S_length) # Numerical index of genes whose means are permitted to vary
  r=2 # Value of r (for negative binomial parameter) is fixed
  Zck  = diag(dim(X)[2])>0 # Each cell is assigned to it's own class

  for(iter in 0:200){
    #if(iter!=0){
    #  Zck = split_clusters
    #}

    # EM steps
    expOut = expectation_step(cl,X,Zck,r,S)
    dim(expOut$Zck)
    dim(expOut$log_Pr_ck_matrix)
    Zck = expOut$Zck
    log_Pr_ck_matrix = expOut$log_Pr_ck_matrix
    S   = maximisation_step(cl,X,Zck,r,S,S_length)

    # Calculate BIC penalty
    BIC_penalty = (S_length * log(apply(Zck,2,sum)))/2

    prune_out = prune_clusters(Zck,log_Pr_ck_matrix,S,BIC_penalty)
    if(prune_out$dropped!=-1){
      Zck = prune_out$Zck
      print(sprintf("------OUTPUT: penalty: %s, diff: %s",prune_out$BIC_penalty,prune_out$diffLogLik))
    }
  }


}
