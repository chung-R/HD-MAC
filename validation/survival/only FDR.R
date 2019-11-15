# data
library(cp4p)
library(data.table)
setwd("C:/Users/Judy/Desktop/Meeting/RShiny package/RShiny")
source("c_index.R")
cox = as.data.frame(fread("Ova_TCGA_OS_clinical_muta_cleaned_313_13_670.csv",header = T,sep = ","))
library(survival)
# 參數設置
user.data = cox
time.col = 9
event.col = 10
gene.type = list(gene.var = colnames(cox)[14:683])
fdr.threshold = 0.05
penalty.alpha = 1
nfold = 1

#----------------------------------------------
# validation
set.seed(1)
t       = sample(1:nfold,size = nrow(user.data), replace = T)
gene.col             = gene.type$gene.var
user.data[,gene.col] = apply(user.data[,gene.col], 2, as.numeric)
time                 = as.numeric(user.data[,time.col])
delta                = user.data[,event.col]
adaptive.lasso       = (penalty.alpha == 2)
gc()

#### Training/Testing ####
train  = !(t==nfold)
test   =  (t==nfold)
if(nfold == 1){train  = (t==nfold);test   =  (t==nfold)} # nfold = 1

#### Delete the 0 sum gene covariate column ####
tmpsum = apply(user.data[ train, gene.col ], 2, sum)
gene.col.train = gene.col[tmpsum != 0]

gene.p.val = sapply(gene.col.train , function(i){
  # fit   = coxph(Surv(time, delta)~., data = as.data.frame(cbind(time = time, delta = delta, user.data[,i])), subset = train)
  fit  = survdiff(Surv(time, delta) ~ ., data = as.data.frame(cbind(time = time, delta = delta, user.data[,i])), subset = train)
  return(1-pchisq(fit$chisq,1)) # summary(fit)$coefficients[5]
})

f <- adjust.p(gene.p.val, pi0.method = 1)
f <- f$adjp[,2]

# result
gene.col.train[ which( f <= fdr.threshold ) ]
