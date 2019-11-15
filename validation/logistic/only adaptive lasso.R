## setting
library(selectiveInference)
library(glmnet)
library(data.table)
setwd("C:/Users/Judy/Desktop/Meeting/RShiny package/RShiny")
lr = as.data.frame(fread("Blad_SUB_8clinical_4937mut_8024mRNA.csv",header = T,sep = ","))
fdr.genes.col = read.csv("C:/Users/Judy/Desktop/Meeting/RShiny package/logistic/fdr_0.05_result.csv")[,-1]
fdr.genes.col = as.character(fdr.genes.col)
user.data = lr
response.col = 3
gene.type = list(gene.var = colnames(lr)[4948:12971]) 
nfold = 1

response   = as.factor(user.data[,response.col]) 
x <- model.matrix( response ~ . , data.frame( response = as.numeric(response == 1), user.data[,fdr.genes.col]) )[,-1]
y <- as.numeric(response == 1)

set.seed(2)
t       = sample(1:nfold,size = nrow(user.data), replace = T)
train  = (t==1);test   =  (t==1) # fold=1
x.train = x[train,]
x.test  = x[test ,]

y.train = y[train]
y.test  = y[test ]
#-------------------------------------------------
library(ncvreg)
cv.out <- cv.ncvreg((x.train), factor(y.train), penalty="lasso", nfolds = 5, family = "binomial", seed = 6)
lambda_min   <- cv.out$lambda.min
if(sum(as.numeric(coef(cv.out, s = lambda_min))==0) == length(fdr.genes.col)){
  lambda_min <- cv.out$lambda[as.numeric(which(cv.out$nzero != 0)[1])]
}
cv.est.coef       = coef(cv.out, s=cv.out$lambda.min) 
best_coef         = as.numeric(cv.est.coef)
geneSel           = rownames(as.matrix(cv.est.coef))[which(!( best_coef == 0))][-1]

#### Weight Change ####
w          = 1 / (abs(best_coef)[which(best_coef != 0)][-1])
w[w == Inf] = 999999999

#### Perform adaptive LASSO ####
fit1_cv                   = cv.glmnet(x = x.train[,colnames(x.train) %in% geneSel], y = y.train, family = "binomial", nfold = 5, alpha = 1, penalty.factor = w)
lambda_min                = fit1_cv$lambda.min
cv.adap.lasso.est.coef    = coef(fit1_cv, s = fit1_cv$lambda.min)
best_coef_adap_la         = as.numeric(cv.adap.lasso.est.coef)
adap_la_geneSel           = rownames(as.matrix(cv.adap.lasso.est.coef))[which(!( best_coef_adap_la == 0))][-1]

#### Standard error and p-value for estimated coefficients ####
  adap.lasso.est.coef       = as.numeric(coef(fit1_cv, s = fit1_cv$lambda.min))
  
  cv.adap.lasso.est.coef    = coef(fit1_cv, s = fit1_cv$lambda.min)
  best_coef_adap_la         = as.numeric(cv.adap.lasso.est.coef)
  adap.l.i                  = fixedLassoInf(x = x[,colnames(x) %in% geneSel], y = y, beta = cv.adap.lasso.est.coef, lambda = fit1_cv$lambda.min, family = "binomial")
  best_coef_adap_la_p_value = adap.l.i$pv
  
  la_geneSel          = row.names(cv.adap.lasso.est.coef)[which(best_coef_adap_la != 0)][-1]
  d = cbind(gene_list = la_geneSel,
            estimated_coefficient = best_coef_adap_la[which(best_coef_adap_la != 0)][-1], 
            p_value = best_coef_adap_la_p_value)


#### 0.01-0.99 threshold find cut point by Youden index ####
lasso_prob        = predict( fit1_cv , newx = x.train[,colnames(x.train) %in% geneSel] , s=lambda_min , type="response")
training_sen = rep(0, 9999)
training_spe = rep(0, 9999)
for(i in 1:9999){
  training_response_by_threshold = rep(0, length(lasso_prob))
  over_threshold_index = which(lasso_prob > 0.0001*i)
  training_response_by_threshold[over_threshold_index] = 1
  
  ## Confusion Table
  # training_table  = table(y.train, training_response_by_threshold)
  training_sen_index = sum((training_response_by_threshold == 1)*(y.train == 1))
  training_spe_index = sum((training_response_by_threshold == 0)*(y.train == 0))
  training_sen[i] = training_sen_index/sum(y.train == 1)
  training_spe[i] = training_spe_index/sum(y.train == 0)
}
#### Youden index ####
Youden_index = rep(0, 9999)
for(i in 1:9999){
  Youden_index[i] = training_sen[i] + training_spe[i] - 1
}
#### Cut Point By Youden index ####
cutpt_index = which(Youden_index == max(Youden_index))
cutpt       = 0.0001*median(cutpt_index) # max(cutpt_index)

#### Testing data ####
pred.prob         = predict(fit1_cv, newx = x.test[,colnames(x.test) %in% geneSel], s=lambda_min, type="response") 

#### use training cut point to find sen spe accuracy ####
testing_response_by_cutpt = rep(0, length(pred.prob))
over_cutpt_index = which(pred.prob > cutpt)
testing_response_by_cutpt[over_cutpt_index] = 1

testing_sen_result = sum((testing_response_by_cutpt == 1)*(y.test == 1))/sum(y.test == 1)
testing_spe_result = sum((testing_response_by_cutpt == 0)*(y.test == 0))/sum(y.test == 0)

#### 0.01-0.99 threshold find cut point by Youden index ####
testing_sen = rep(0, 9999)
testing_spe = rep(0, 9999)
for(i in 1:9999){
  testing_response_by_threshold = rep(0, length(pred.prob))
  over_threshold_index = which(pred.prob > 0.0001*i)
  testing_response_by_threshold[over_threshold_index] = 1
  
  ## Confusion Table
  testing_sen_index = sum((testing_response_by_threshold == 1)*(y.test == 1))
  testing_spe_index = sum((testing_response_by_threshold == 0)*(y.test == 0))
  testing_sen[i] = testing_sen_index/sum(y.test == 1)
  testing_spe[i] = testing_spe_index/sum(y.test == 0)
}

#### Calculating testing AUC from ROC curve ####
auc_x  = 1-testing_spe
auc_y  = testing_sen
auc_df = as.data.frame(cbind(auc_x,auc_y))
order_df = auc_df[order(auc_df$auc_x,auc_df$auc_y),]

upper_vector = order_df$auc_y[-dim(order_df)[1]]
lower_vector = order_df$auc_y[-1]
upper_lower_sum = upper_vector + lower_vector
testing_auc  = sum(upper_lower_sum*diff(order_df$auc_x)/2)
rocobj       = list(auc = testing_auc)

#### Print result ####
la1 <- list( cp = cutpt,
             testing_pre0_res0 = sum((testing_response_by_cutpt == 0)*(y.test == 0)),
             testing_pre0_res1 = sum((testing_response_by_cutpt == 0)*(y.test == 1)),
             testing_pre1_res0 = sum((testing_response_by_cutpt == 1)*(y.test == 0)),
             testing_pre1_res1 = sum((testing_response_by_cutpt == 1)*(y.test == 1)),
             test.sen      = testing_sen_result,
             test.spe      = testing_spe_result,
             test.auc      = as.numeric(rocobj$auc),
             selected.gene = (as.matrix(cv.adap.lasso.est.coef))[which(as.matrix(cv.adap.lasso.est.coef) != 0),],
             coef.and.p    = as.data.frame(d)
)

# result
la1$selected.gene
