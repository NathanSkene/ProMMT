prune_clusters <- function(X,Zck,r,S,BIC_penalty){

  log_Pr_ck_matrix = get_log_Pr_ck_matrix(cl,X,Zck,r,S)

  #for(c in 1:dim(Zck)[2]){
  for(c in sample(1:dim(Zck)[2],dim(Zck)[2])){
    # Find which clusters each cell would be assigned to if not the current one
    alt_log_Pr_ck_matrix = log_Pr_ck_matrix
    alt_log_Pr_ck_matrix[Zck[,c],c] = -Inf
    Kck = get_best_clust_given_logPr_ck(alt_log_Pr_ck_matrix)>0
    alt_log_Pr_ck_matrix = get_log_Pr_ck_matrix(cl,X,Kck,r,S) #< Regenerate the whole thing to account for BIC penalty
    # Add the BIC penalty to alternative logPr
    #altBIC_penalty = (S_length * log(apply(Kck,2,sum)))/2
    #alt_log_Pr_ck_matrix = t(t(alt_log_Pr_ck_matrix)+altBIC_penalty)

    # For each cell in the cluster
    #for(cc in which(Zck[,c])){
    #  curLogLik = log_Pr_ck_matrix[cc,Zck[cc,]]
    #  altLogLik = log_Pr_ck_matrix[cc,Kck[cc,]]
    #  diffLogLik = curLogLik - altLogLik
    #}

    # Find current log likelihood
    #curLogLik  = apply(log_Pr_ck_matrix[which(Zck[,c]),],1,max)
    #curLogLik = log_Pr_ck_matrix[which(Zck[,c])]
    curLogLik = log_Pr_ck_matrix[Zck[,c],c]

    # Find what log likelihood would be using the alternative clusters
    #altLogLik = apply(alt_log_Pr_ck_matrix[which(Zck[,c]),],1,max)
    altLogLik  = alt_log_Pr_ck_matrix[which(Zck[,c]),][Kck[which(Zck[,c]),]]

    # Difference
    diffLogLik = sum((altLogLik-curLogLik))

    print(sprintf("Looking at removing %s... Penalty: %s, Diff: %s",c,BIC_penalty[c],diffLogLik))

    # Is the loss of likelihood from merging smaller than the BIC penalty?
    #if(altLogLik>curLogLik){makeChange=TRUE
    if(diffLogLik>=0){makeChange=TRUE  #> Positive diff indidicate the change improves log likelihood... so always make that change
    }#else{makeChange=FALSE}
    else{
       if(diffLogLik < -BIC_penalty[c]){
         makeChange=FALSE
       }else{
         makeChange=TRUE
       }
    }
    #if(diffLogLik<BIC_penalty[c]){
    if(makeChange){
      print("Removing")
      out = list()
      out$dropped = c
      out$BIC_penalty = BIC_penalty[c]
      out$num_cells_changed = sum(Zck[,c])
      out$diffLogLik = diffLogLik
      #out$Zck = Kck[,setdiff(1:dim(Kck)[2],c)]
      out$Zck = Kck[,colSums(Kck)>0]
      return(out)
    }
  }
  # If no clusters are altered that's it
  out = list()
  out$dropped = -1
  return(out)
}
