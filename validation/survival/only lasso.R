# setting
library(data.table)
setwd("C:/Users/Judy/Desktop/Meeting/RShiny package/RShiny")
source("c_index.R")
cox = as.data.frame(fread("Ova_TCGA_OS_clinical_muta_cleaned_313_13_670.csv",header = T,sep = ","))
library(survival)

user.data = cox
time.col = 9
event.col = 10
gene.type = list(gene.var = colnames(cox)[14:683])
penalty.alpha = 1
nfold = 1
fdr.genes.col = c("ZSWIM8", "PABPC3")

set.seed(1)
t       = sample(1:nfold,size = nrow(user.data), replace = T)
gene.col             = gene.type$gene.var
user.data[,gene.col] = apply(user.data[,gene.col], 2, as.numeric)
time                 = as.numeric(user.data[,time.col])
delta                = user.data[,event.col]
x <- data.matrix(user.data[, fdr.genes.col])
y <- Surv(time,delta)

train  = !(t==nfold)
test   =  (t==nfold)
if(nfold == 1){train  = (t==nfold);test   =  (t==nfold)} # nfold = 1

x.train = x[train,]
x.test  = x[test ,]
y.train = with(data.frame(time,delta)[train,], Surv(as.numeric(time),delta))
y.test  = with(data.frame(time,delta)[test, ], Surv(as.numeric(time),delta))

#-----------------------------------------------------
# validation
library(ncvreg)
set.seed(4)
cv.out1 <- cv.ncvsurv(x.train, y.train, penalty="lasso", nfolds = 5)

lambda_min1   <- cv.out1$lambda.min
if(sum(as.numeric(coef(cv.out1, s = lambda_min1))==0) == length(fdr.genes.col)){
  nz <- unlist(lapply(c(1:length(cv.out1$lambda)), function(x) sum(coef(cv.out1, s = cv.out1$lambda[x])!=0))) # cv.out$nzero
  lambda_min1 <- cv.out1$lambda[as.numeric(which(nz != 0)[1])]
}

lasso.est.coef1    = as.numeric(coef(cv.out1, s = lambda_min1))
cv.est.coef1       = coef(cv.out1, s = lambda_min1/dim(user.data)[1])
best_coef1         = as.numeric(cv.est.coef1)
li1 = fixedLassoInf(x, time, beta = best_coef1, lambda = lambda_min1, status = delta, family = "cox")

best_coef_p_value1 = li1$pv

cv.lasso.est.coef1   = coef(cv.out1, s = lambda_min1)
la_geneSel1          = row.names(cv.est.coef1)[which(best_coef1 != 0)]
d1 = cbind(gene_list = la_geneSel1, 
           estimated_coefficient = lasso.est.coef1[which(lasso.est.coef1 != 0)], 
           p_value = best_coef_p_value1)

u51          =  data.frame(Death = delta, Death_surtime = time)[test,]
cut_hazard11 =  predict(cv.out1 , X = x.test, s=lambda_min1)

la1 <- list( c_index       = c.index(u5 = u51, cut_hazard1 = cut_hazard11),
             selected.gene = ((as.matrix(coef(cv.out1, s = lambda_min1)))[which( (as.numeric(coef(cv.out1, s = lambda_min1))) != 0),] ),
             coef.and.p    = as.data.frame(d1) )

# result
la1$selected.gene
