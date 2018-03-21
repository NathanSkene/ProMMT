# Expectation step
expectation_step <- function(cl,X,Zck,r,S){

  log_Pr_ck_matrix = get_log_Pr_ck_matrix(cl,X,Zck,r,S)

  Zck_new = get_best_clust_given_logPr_ck(log_Pr_ck_matrix)

  # If a cluster now has no cells, then drop that cluster
  keepCols = !colSums(Zck_new)==0
  Zck_new = Zck_new[,keepCols]

  # Force Zck to be logical
  Zck_new = Zck_new>0

  return(list(Zck=Zck_new))
}
