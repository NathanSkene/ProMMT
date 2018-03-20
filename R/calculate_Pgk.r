# Calculate negative binomial p parameters
calculate_Pgk <- function(Ugk,r){
  Pgk = Ugk/(r+Ugk)
  return(Pgk)
}
