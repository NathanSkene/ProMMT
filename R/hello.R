library(parallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores-1)


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
X=exp
#res = promtt(exp,cl)

S_length = 150 # How many genes shoul be permitted to have different means?
S = sample(1:dim(X)[1],S_length) # Numerical index of genes whose means are permitted to vary
r=2 # Value of r (for negative binomial parameter) is fixed
Zck  = diag(dim(X)[2])>0 # Each cell is assigned to it's own class
