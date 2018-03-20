get.deltaYgt <- function(g,Xc,r){
  regularised_mean <- function(x){   out = (10^-4 + sum(x)) / (1+length(x))   }
  Xcg = Xc[g,]
  thresholds = sort(unique(Xc[g,]))[-1]
  for(t in thresholds){
    Xcg_above = Xcg[Xcg>=t]
    Xcg_below = Xcg[Xcg<t]
    meanAboveThresh = regularised_mean(Xcg_above)
    meanBelowThresh = regularised_mean(Xcg_below)
    Pg        = regularised_mean(Xcg)/(r+regularised_mean(Xcg))
    Pgt_above = meanAboveThresh/(r+meanAboveThresh)
    Pgt_below = meanBelowThresh/(r+meanBelowThresh)
    deltaYgt_aboveTerm = sum(Xcg_above*(log(Pgt_above)-log(Pg))  +  r*(log(1-Pgt_above)-log(1-Pg)))
    deltaYgt_belowTerm = sum(Xcg_below*(log(Pgt_below)-log(Pg))  +  r*(log(1-Pgt_below)-log(1-Pg)))
    deltaYgt = deltaYgt_aboveTerm + deltaYgt_belowTerm
    tmp = data.frame(g=g,t=t,deltaYgt=deltaYgt)
    if(t==thresholds[1]){out=tmp}else{out=rbind(out,tmp)}
  }
  if(length(thresholds)>0){   return(out)   }
}
