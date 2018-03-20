calculate_priors <- function(Zck){
  # Priors are set as the number of cells in each cluster
  return(apply(Zck,2,sum)/dim(Zck)[1])
}
