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
library(penalized)
aa <- optL2(y.train, x.train, lambda1 = 0, fold = 5, model = "cox")
aa2 <- penalized(y.train, x.test, lambda2=aa$lambda,
                 lambda1=0, model = "cox")

u5          =  data.frame(Death = delta, Death_surtime = time)[test,]
cut_hazard1 = matrix(as.matrix(predict(aa2,penalized=x.test))[,198], ncol = 1)

#### print result ####

la1 <- list( c_index = c.index(u5 = u5, 
                               cut_hazard1 = cut_hazard1), 
             selected.gene = ((as.matrix(coefficients(aa2), ncol = 1))[which( (as.numeric(coefficients(aa2))) != 0),] ))

# result
la1$selected.gene
