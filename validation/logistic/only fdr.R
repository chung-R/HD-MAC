## FDR
# setting
setwd("C:/Users/Judy/Desktop/Meeting/RShiny package/RShiny")
lr = as.data.frame(fread("Blad_SUB_8clinical_4937mut_8024mRNA.csv",header = T,sep = ","))
library(cp4p)
user.data = lr
response.col = 3
gene.type = list(gene.var = colnames(lr)[4948:12971]) 
nfold = 1
fdr.threshold  = 0.05

set.seed(2)
t       = sample(1:nfold,size = nrow(user.data), replace = T)
gene.col = gene.type$gene.var
response = as.factor(user.data[,response.col]) 
train  = (t==1);test   =  (t==1) # train and test
#### Delete the 0 sum gene covariate column ####
tmpsum = apply(user.data[ train, gene.col ], 2, sum)
gene.col.train = gene.col[tmpsum != 0]

#### univariate logistic regression fdr ####
tmp.p = rep(NA, length(gene.col.train))
tmp.p = sapply( gene.col.train, function(i){ summary( glm( response ~ . , family=binomial, data = data.frame(response = as.numeric(response == 1), as.numeric(user.data[,i])), subset = train ) )$coefficients[2,4] })
f <- adjust.p(tmp.p, pi0.method = 1)
f <- f$adjp[,2]

# result
gene.col.train[ which( f <= fdr.threshold ) ]
