get_best_clust_given_logPr_ck <- function(log_Pr_ck_matrix){
  # Now select new values of Zck (cells-->clusters based on the maximum value of Pck in each row)
  Zck_new <- t(apply(log_Pr_ck_matrix, 1, function(z){
    1 * (z == max(z))
  }))

  # Check only one cluster per row has 1
  for(ro in which(rowSums(Zck_new)>1)){
    whichDup = duplicated(Zck_new[ro,]) & Zck_new[ro,]==1
    Zck_new[ro,whichDup] = rep(0,length(sum(whichDup)))
  }
  return(Zck_new)
}
