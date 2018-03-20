

# How is S defined? S, meanK and pi are fit my maximum likelihood

# Calculate the probability of gene expression vector for a single cell: formula 1
#probability_of_x <- function(X,Ugk){
#}

# Calculate negative binomial p parameters
calculate_Pgk <- function(Ugk,r){
  Pgk = Ugk/(r+Ugk)
  return(Pgk)
}

# Calculate the probability of gene expression vector for a single cell belonging to k: formula 2
probability_of_x_given_k <- function(X,k,S){

}

# Probability of read value x given negative binomial with r=2 and p=?
probability_of_x_given_r_p <-function(x,r,p){
  # Function 3
}

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

# Maximisation step
maximisation_step <- function(cl,X,Zck,r,S,S_length){
  # S and u are optimised given current values of z_c_k
  Ugk  = calculate_cluster_mean(X,Zck,cl)
  Ug0  = apply(X,1,mean)
  Pg0  = Ug0/(r+Ug0)
  Pgk  = calculate_Pgk(Ugk,r)

  # Calculate for each gene, Yg (the gain in log likelihood when the distribution of gene g is allowed to vary between classes)
  calc_Yg <- function(g,Zck,X,Pgk,Pg0,r){
    val = 0
    for(c in 1:dim(Zck)[2]){
      val = val + X[g,c]*(log(Pgk[g,c] - log(Pg0[g]))) + r*(log(1-Pgk[g,c]) - log(1-Pg0[g]))
    }
    return(val)
  }
  Yg = unlist(parLapply(cl,as.list(1:dim(X)[1]),FUN=calc_Yg,Zck,X,Pgk,Pg0,r))

  # Select the new set of S (those with biggest values of Yg)
  top_Yg_genes = names(sort(Yg,decreasing=TRUE)[1:S_length])
  S_new = which(rownames(X) %in% top_Yg_genes)
  return(S_new)
}

split_clusters <- function(){

}

prune_clusters <- function(Zck,log_Pr_ck_matrix,S,BIC_penalty){
  #for(c in 1:dim(Zck)[2]){
  for(c in sample(1:dim(Zck)[2],dim(Zck)[2])){
    # Find which clusters each cell would be assigned to if not the current one
    alt_log_Pr_ck_matrix = log_Pr_ck_matrix
    alt_log_Pr_ck_matrix[Zck[,c],c] = -Inf
    Kck = get_best_clust_given_logPr_ck(alt_log_Pr_ck_matrix)>0

    # For each cell in the cluster
    for(cc in which(Zck[,c])){
      curLogLik = log_Pr_ck_matrix[cc,Zck[cc,]]
      altLogLik = log_Pr_ck_matrix[cc,Kck[cc,]]
      diffLogLik = curLogLik - altLogLik
    }

    # Find current log likelihood
    curLogLik = log_Pr_ck_matrix[which(Zck[,c])]

    # Find what log likelihood would be using the alternative clusters
    altLogLik = log_Pr_ck_matrix[which(Zck[,c]),][which(Kck[which(Zck[,c]),])]

    # Difference
    diffLogLik = sum((altLogLik-curLogLik))

    print(sprintf("Looking at removing %s... Penalty: %s, Diff: %s",c,BIC_penalty[c],diffLogLik))

    # Is the loss of likelihood from merging smaller than the BIC penalty?
    if(altLogLik>curLogLik){makeChange=TRUE
    }else{
      if(diffLogLik < -BIC_penalty[c]){
        makeChange=TRUE
      }
      makeChange=FALSE
    }
    #if(diffLogLik<BIC_penalty[c]){
    if(makeChange){
      print("Removing")
      out = list()
      out$dropped = c
      out$BIC_penalty = BIC_penalty[c]
      out$num_cells_changed = sum(Zck[,c])
      out$diffLogLik = diffLogLik
      out$Zck = Kck[,setdiff(1:dim(Kck)[2],c)]
      return(out)
    }
  }
  # If no clusters are altered that's it
  out = list()
  out$dropped = -1
  return(out)
}

calculate_cluster_mean <- function(X,Zck,cl){
  # Error check
  if(class(Zck[,1])!="logical"){stop("Zck must be logical, not numeric")}

  get.mean.geneexp <- function(whichCells,X){
      regularised_mean <- function(x){
        out = (10^-4 + sum(x)) / (1+length(x))
      }
      print(which(whichCells))
      #print(X[,whichCells,drop=FALSE])
      return(apply(X[,whichCells,drop=FALSE],1,regularised_mean))
  }
  Ugk = parApply(cl,Zck,2,get.mean.geneexp,X)

  return(Ugk)
}

calculate_priors <- function(Zck){
  # Priors are set as the number of cells in each cluster
  return(apply(Zck,2,sum)/dim(Zck)[1])
}

# Run the ProMTT algorithm
ProMTT <- function(X,cl){
  # Initial setup
  S_length = 150 # How many genes shoul be permitted to have different means?
  S = sample(1:dim(X)[1],S_length) # Numerical index of genes whose means are permitted to vary
  r=2 # Value of r (for negative binomial parameter) is fixed
  Zck  = diag(dim(X)[2])>0 # Each cell is assigned to it's own class

  for(iter in 0:100){
    #if(iter!=0){
    #  Zck = split_clusters
    #}

    # EM steps
    expOut = expectation_step(cl,X,Zck,r,S)
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

library(parallel)
## Use option cl.cores to choose an appropriate cluster size.
no_cores <- detectCores()
cl <- makeCluster(no_cores)

# Load dummy data
library(EWCE)
data("cortex_mrna")
annot = cortex_mrna$annot[cortex_mrna$annot$level2class %in% c("Int10","Mgl1"),]
exp   = cortex_mrna$exp[,annot$cell_id]
library(limma)
design <- model.matrix(~as.character(annot$level2class))
fit <- lmFit(exp, design)
fit <- eBayes(fit)
tt  <- topTable(fit, coef = "as.character(annot$level2class)Mgl1", adjust = "BH",number=10000000)
ttU = tt[order(tt$t,decreasing=TRUE),]
ttD = tt[order(tt$t,decreasing=FALSE),]
up500=rownames(ttU[1:500,])
down500=rownames(ttD[1:500,])
keepGenes = c(up500,down500)
exp = exp[keepGenes,]
res = promtt(exp,cl)
