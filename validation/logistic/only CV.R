## setting
library(selectiveInference)
library(glmnet)
library(cvTools)
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

#t11--> cv.glmnet 用cvTools::cvFolds取代sample
t11 <- function (x, y, weights, offset = NULL, lambda = NULL, type.measure = c("mse", 
                                                                               "deviance", "class", "auc", "mae"), 
                 nfolds = 10, foldid, alignment = c("lambda", "fraction"), 
                 grouped = TRUE, keep = FALSE, parallel = FALSE, ...) 
{
  if (missing(type.measure)) 
    type.measure = "default"
  else type.measure = match.arg(type.measure)
  alignment = match.arg(alignment)
  if (!is.null(lambda) && length(lambda) < 2) 
    stop("Need more than one value of lambda for cv.glmnet")
  if (!is.null(lambda) && alignment == "fraction") {
    warning("fraction of path alignment not available if lambda given as argument; switched to alignment=`lambda`")
    alignment = "lambda"
  }
  N = nrow(x)
  if (missing(weights)) 
    weights = rep(1, N)
  else weights = as.double(weights)
  y = drop(y)
  glmnet.call = match.call(expand.dots = TRUE)
  which = match(c("type.measure", "nfolds", "foldid", 
                  "grouped", "keep"), names(glmnet.call), F)
  if (any(which)) 
    glmnet.call = glmnet.call[-which]
  glmnet.call[[1]] = as.name("glmnet")
  glmnet.object = glmnet(x, y, weights = weights, offset = offset, 
                         lambda = lambda, ...)
  glmnet.object$call = glmnet.call
  subclass = class(glmnet.object)[[1]]
  type.measure = cvtype(type.measure, subclass)
  is.offset = glmnet.object$offset
  if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
    nz = predict(glmnet.object, type = "nonzero")
    nz = sapply(nz, function(x) sapply(x, length))
    nz = ceiling(apply(nz, 1, median))
  }
  else nz = sapply(predict(glmnet.object, type = "nonzero"), 
                   length)
  if (missing(foldid)) {
    t1 <- cvFolds(n = N, K = nfolds)
    foldid = t1$which[order(t1$subsets)]}
  else nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist = as.list(seq(nfolds))
  if (parallel) {
    outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar% 
    {
      which = foldid == i
      if (length(dim(y)) > 1) 
        y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset) 
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      glmnet(x[!which, , drop = FALSE], y_sub, lambda = lambda, 
             offset = offset_sub, weights = weights[!which], 
             ...)
    }
  }
  else {
    for (i in seq(nfolds)) {
      which = foldid == i
      if (is.matrix(y)) 
        y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset) 
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      outlist[[i]] = glmnet(x[!which, , drop = FALSE], 
                            y_sub, lambda = lambda, offset = offset_sub, 
                            weights = weights[!which], ...)
    }
  }
  fun = paste("cv", subclass, sep = ".")
  lambda = glmnet.object$lambda
  cvstuff = do.call(fun, list(outlist, lambda, x, y, weights, 
                              offset, foldid, type.measure, grouped, keep, alignment))
  cvm = cvstuff$cvm
  cvsd = cvstuff$cvsd
  nas = is.na(cvsd)
  if (any(nas)) {
    lambda = lambda[!nas]
    cvm = cvm[!nas]
    cvsd = cvsd[!nas]
    nz = nz[!nas]
  }
  cvname = names(cvstuff$type.measure)
  names(cvname) = cvstuff$type.measure
  out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + 
               cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, glmnet.fit = glmnet.object)
  if (keep) 
    out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
  lamin = if (cvname == "AUC") 
    getmin(lambda, -cvm, cvsd)
  else getmin(lambda, cvm, cvsd)
  obj = c(out, as.list(lamin))
  class(obj) = "cv.glmnet"
  obj
}



cv.out1       <- t11(x.train, y.train, alpha = penalty.alpha, family="binomial", type.measure = "mse", nfolds = 5 )
#-----------------------------------------------
lambda_min   <- cv.out1$lambda.min
if(sum(as.numeric(coef(cv.out1, s = lambda_min))==0) == length(fdr.genes.col)){
  lambda_min <- cv.out1$lambda[as.numeric(which(cv.out1$nzero != 0)[1])]
}  
#### Standard error and p-value for estimated coefficients ####

#### training data: find cut point for testing data ####
lasso_prob1 = predict(cv.out1 , newx = x.train , s = lambda_min , type = "response")

# result
lasso_prob1
