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
library(glmnet)
set.seed(4)
cv.out1 <- cv.ncvsurv(x.train, y.train, penalty="lasso", nfolds = 5)

lambda_min1   <- cv.out1$lambda.min
if(sum(as.numeric(coef(cv.out1, s = lambda_min1))==0) == length(fdr.genes.col)){
  nz <- unlist(lapply(c(1:length(cv.out1$lambda)), function(x) sum(coef(cv.out1, s = cv.out1$lambda[x])!=0))) # cv.out$nzero
  lambda_min1 <- cv.out1$lambda[as.numeric(which(nz != 0)[1])]
}



#### Adaptive Lasso : Weight Change ####
cv.est.coef       = coef(cv.out1, s =  lambda_min  )
best_coef         = as.numeric(coef(cv.out1, s =  lambda_min  ))
geneSel           = rownames(as.matrix(cv.est.coef))[which(!( best_coef == 0))]
w.vector          = 1 / (abs(best_coef)[which(!( best_coef == 0))])
w.vector[w.vector == Inf] = 999999999

#### Perform adaptive LASSO ####
fit1_cv           = cv.glmnet(x = x.train[,colnames(x.train) %in% geneSel], y = y.train, family = "cox",nfold = 5, alpha = 1, penalty.factor = w.vector, type.measure = "deviance")
lambda_min_ada    = fit1_cv$lambda.min
cv.adap.lasso.est.coef    = coef(fit1_cv, s = fit1_cv$lambda.min)
best_coef_adap_la         = as.numeric(cv.adap.lasso.est.coef)
adap_la_geneSel           = rownames(cv.adap.lasso.est.coef)

#### Standard error and p-value for estimated coefficients ####
if( nfold == 1 && penalty.alpha == 1 ){
  adaplasso.est.coef        = as.numeric(coef(fit1_cv, s = lambda_min_ada))
  
  cv.adap.lasso.est.coef    = coef(fit1_cv, s = lambda_min_ada/dim(user.data)[1])
  best_coef_adap_la         = as.numeric(cv.adap.lasso.est.coef)
  adap.l.i                  = fixedLassoInf(x = x[,colnames(x) %in% geneSel], time, beta = best_coef_adap_la, lambda = lambda_min_ada, status = delta, family = "cox")
  # dd = x.train[,colnames(x.train) %in% geneSel]
  # index = which(duplicated(x.train[,colnames(x.train) %in% geneSel], MARGIN = 2))
  # adap.l.i = fixedLassoInf(dd[,-index], time, beta = best_coef_adap_la[-index], lambda = lambda_min, status = delta, family = "cox")
  # which(duplicated(x, MARGIN = 2))
  best_coef_adap_la_p_value = adap.l.i$pv
  
  la_geneSel        = row.names(cv.adap.lasso.est.coef)[which(best_coef_adap_la != 0)]
  d = cbind(gene_list = la_geneSel,
            estimated_coefficient = adaplasso.est.coef[which(adaplasso.est.coef != 0)], 
            p_value = best_coef_adap_la_p_value)
}

lasso_prob        = predict( fit1_cv , newx = x.train[,colnames(x.train) %in% geneSel], s=lambda_min_ada , type="response")

#### Testing data ####
u5                = data.frame(Death = delta, Death_surtime = time)[test,]
cut_hazard1       = predict(fit1_cv, newx = x.test[,colnames(x.test) %in% geneSel], s=lambda_min_ada)

#### print result ####
la <- list( c_index       = c.index(u5 = u5,cut_hazard1 = cut_hazard1),
            selected.gene =( (as.matrix(coef(fit1_cv, s = lambda_min_ada))) )[which((as.numeric(coef(fit1_cv, s = lambda_min_ada))) != 0), ],
            coef.and.p    = as.data.frame(d) )
   
la$selected.gene

