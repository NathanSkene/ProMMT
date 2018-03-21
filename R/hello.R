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
#res = promtt(exp,cl)
