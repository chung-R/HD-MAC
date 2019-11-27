#######################################################################
# Cox PH Model 
#######################################################################
library(survival)
library(glmnet)
library(data.table)
library(selectiveInference)

#################################################

cox.gene.penalty.cv.result = function(user.data, time.col, event.col, 
                                      gene.type, padj.method   = "fdr", fdr.threshold = 0.05,
                                      penalty.alpha = 1, nfold = 1, force.add.cov = NULL ){
  set.seed(1)
  t       = sample(1:nfold,size = nrow(user.data), replace = T)
  gene.col             = gene.type$gene.var
  user.data[,gene.col] = apply(user.data[,gene.col], 2, as.numeric)
  time                 = as.numeric(user.data[,time.col])
  delta                = user.data[,event.col]
  adaptive.lasso       = (penalty.alpha == 2)
  gc()
  
  #### Training/Testing ####
  cox.gene.penalty.cv.list = sapply(1:nfold, function(fold){
    train  = !(t==fold)
    test   =  (t==fold)
    if(nfold == 1){train  = (t==fold);test   =  (t==fold)} # fold = 1
    
    #### Delete the 0 sum gene covariate column ####
    tmpsum = apply(user.data[ train, gene.col ], 2, sum)
    gene.col.train = gene.col[tmpsum != 0]
    
    #### univariate logistic regression fdr ####
    
    if( padj.method == "fdr" ){ 
      gene.p.val = sapply(gene.col.train , function(i){
        # fit   = coxph(Surv(time, delta)~., data = as.data.frame(cbind(time = time, delta = delta, user.data[,i])), subset = train)
        fit  = survdiff(Surv(time, delta) ~ ., data = as.data.frame(cbind(time = time, delta = delta, user.data[,i])), subset = train)
        return(1-pchisq(fit$chisq,1)) # summary(fit)$coefficients[5]
      })
      tmp.p = p.adjust(gene.p.val, padj.method)
      fdr.genes.col = gene.col.train[ which( tmp.p <= fdr.threshold ) ]
      if (length(fdr.genes.col) == 0) stop("small fdr.threshold --> no genes selected")}
    else{
        fdr.genes.col = gene.col.train
      }
    
    #### Regression part ####
    x <- data.matrix(user.data[, fdr.genes.col])
    y <- Surv(time,delta)
    
    x.train = x[train,]
    x.test  = x[test ,]
    y.train = with(data.frame(time,delta)[train,], Surv(as.numeric(time),delta))
    y.test  = with(data.frame(time,delta)[test, ], Surv(as.numeric(time),delta))
    
    if( adaptive.lasso ){penalty.alpha = 1}
    
    set.seed(4)
    cv.out       <- cv.glmnet(x.train, y.train, alpha = penalty.alpha, family="cox", type.measure = "deviance", nfolds = 5 )
    gc()
    lambda_min   <- cv.out$lambda.min
    if(sum(as.numeric(coef(cv.out, s = lambda_min))==0) == length(fdr.genes.col)){
      lambda_min <- cv.out$lambda[as.numeric(which(cv.out$nzero != 0)[1])]
      }
    
    #### Ridge and Lasso ####
    if( ! adaptive.lasso ){
      
      #### Standard error and p-value for estimated coefficients ####
      if( nfold == 1 && penalty.alpha == 1 ){
        lasso.est.coef    = as.numeric(coef(cv.out, s = lambda_min))
        cv.est.coef       = coef(cv.out, s = lambda_min/dim(user.data)[1])
        best_coef         = as.numeric(cv.est.coef)
        ##li = fixedLassoInf(x, time, beta = best_coef, lambda = lambda_min, status = delta, family = "cox")
        # li = fixedLassoInf(x[,-which(duplicated(x, MARGIN = 2))], time, beta = best_coef[-which(duplicated(x, MARGIN = 2))], lambda = lambda_min, status = delta, family = "cox")
        # which(duplicated(x, MARGIN = 2))
        ##best_coef_p_value = li$pv
        
        cv.lasso.est.coef   = coef(cv.out, s = lambda_min)
        la_geneSel          = row.names(cv.est.coef)[which(best_coef != 0)]
        d = cbind(gene_list = la_geneSel, 
                  estimated_coefficient = lasso.est.coef[which(lasso.est.coef != 0)] #, 
                  ##p_value = best_coef_p_value
                  )
        }
      
      #### c-index ####
      u5          =  data.frame(Death = delta, Death_surtime = time)[test,]
      cut_hazard1 =  predict(cv.out , newx = x.test, s=lambda_min)
      
      #### print result ####
      if(is.null(force.add.cov)){
        if( nfold == 1 && penalty.alpha == 1 ){
        return( list( c_index       = c.index(u5 = u5, cut_hazard1 = cut_hazard1), 
                      selected.gene = ((as.matrix(coef(cv.out, s = lambda_min)))[which( (as.numeric(coef(cv.out, s = lambda_min))) != 0),] ),
                      coef.and.p    = as.data.frame(d) ))}
        else{return( list( c_index       = c.index(u5 = u5, cut_hazard1 = cut_hazard1), 
                           selected.gene = ((as.matrix(coef(cv.out, s = lambda_min)))[which( (as.numeric(coef(cv.out, s = lambda_min))) != 0),] )))}
      }else{
        #### add clinical covariates ####
        geneSel              = rownames(( (as.matrix(coef(cv.out, s = lambda_min))) ))[which( (as.numeric(coef(cv.out, s = lambda_min))) != 0)]
        cox.data             = data.frame(cbind( time = time[train], delta = delta[train],user.data[train, force.add.cov],x.train[, c(which(colnames(x.train) %in% geneSel)) ]) )
        fit.geneSel.forceCov = coxph( Surv(time, delta) ~ . , data = cox.data)
        
        #### Testing data ####
        cut_hazard1          =  predict(fit.geneSel.forceCov, newdata =  data.frame(cbind( time = time[test], delta = delta[test],user.data[test, force.add.cov],x.test[, c(which(colnames(x.test) %in% geneSel)) ]) ), type="lp")
        u5                   =  data.frame(Death = delta, Death_surtime = time)[test,]
        summary.fit          =  summary(fit.geneSel.forceCov)$coefficients
        rownames(summary.fit)[1:length(force.add.cov)] = colnames(user.data)[force.add.cov]
        
        #### print result ####
        if( nfold == 1 && penalty.alpha == 1 ){
        return(list( c_index = c.index(u5 = u5, cut_hazard1 = cut_hazard1), 
                     selected.gene = (summary.fit[,1])[! is.na(summary.fit[,1])],
                     coef.and.p    = as.data.frame(d) ))}
        else{return(list( c_index = c.index(u5 = u5, cut_hazard1 = cut_hazard1), 
                          selected.gene = (summary.fit[,1])[! is.na(summary.fit[,1])] ))}
      }
    }else{
      #### Adaptive Lasso : Weight Change ####
      cv.est.coef       = coef(cv.out, s =  lambda_min  )
      best_coef         = as.numeric(coef(cv.out, s =  lambda_min  ))
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
        ##adap.l.i                  = fixedLassoInf(x = x[,colnames(x) %in% geneSel], time, beta = best_coef_adap_la, lambda = lambda_min_ada, status = delta, family = "cox")
        # dd = x.train[,colnames(x.train) %in% geneSel]
        # index = which(duplicated(x.train[,colnames(x.train) %in% geneSel], MARGIN = 2))
        # adap.l.i = fixedLassoInf(dd[,-index], time, beta = best_coef_adap_la[-index], lambda = lambda_min, status = delta, family = "cox")
        # which(duplicated(x, MARGIN = 2))
        ##best_coef_adap_la_p_value = adap.l.i$pv
        
        la_geneSel        = row.names(cv.adap.lasso.est.coef)[which(best_coef_adap_la != 0)]
        d = cbind(gene_list = la_geneSel,
                  estimated_coefficient = adaplasso.est.coef[which(adaplasso.est.coef != 0)] #, 
                  ##p_value = best_coef_adap_la_p_value
                  )
        }
      
      lasso_prob        = predict( fit1_cv , newx = x.train[,colnames(x.train) %in% geneSel], s=lambda_min_ada , type="response")
      
      #### Testing data ####
      u5                = data.frame(Death = delta, Death_surtime = time)[test,]
      cut_hazard1       = predict(fit1_cv, newx = x.test[,colnames(x.test) %in% geneSel], s=lambda_min_ada)
      
      #### print result ####
      if(is.null(force.add.cov)){
        if( nfold == 1 && penalty.alpha == 1 ){
          return(list( c_index       = c.index(u5 = u5,cut_hazard1 = cut_hazard1),
                       selected.gene =( (as.matrix(coef(fit1_cv, s = lambda_min_ada))) )[which((as.numeric(coef(fit1_cv, s = lambda_min_ada))) != 0), ],
                       coef.and.p    = as.data.frame(d) ))}
        else{return(list( c_index       = c.index(u5 = u5,cut_hazard1 = cut_hazard1),
                          selected.gene =( (as.matrix(coef(fit1_cv, s = lambda_min_ada))) )[which((as.numeric(coef(fit1_cv, s = lambda_min_ada))) != 0), ]))}
        
      }else{
        #### add clinical covariates ####
        geneSel              = rownames( as.matrix(coef(fit1_cv, s = lambda_min_ada)) )[which( (as.numeric(coef(fit1_cv, s = lambda_min_ada))) != 0)]
        if( length(geneSel) == 0 ) stop("no genes selected through adaptive lasso")
        cox.data             = data.frame(cbind( time = time[train], delta = delta[train], user.data[train, force.add.cov], x.train[, c(which(colnames(x.train) %in% geneSel)) ]) )
        fit.geneSel.forceCov = coxph(Surv(time, delta) ~ ., data = cox.data)
        
        #### Testing data ####
        u5                   = data.frame(Death = delta, Death_surtime = time)[test,]
        cut_hazard1          = predict(fit.geneSel.forceCov, newdata =  data.frame(cbind(time = time[test], 
                                                                                         delta = delta[test], 
                                                                                         user.data[test, force.add.cov],
                                                                                         x.test[, c(which(colnames(x.test) %in% geneSel)) ]) ), 
                                       type="lp")
        
        summary.fit          = summary(fit.geneSel.forceCov)$coefficients
        rownames(summary.fit)[1:length(force.add.cov)] = colnames(user.data)[force.add.cov]
        
        #### print result ####
        if( nfold == 1 && penalty.alpha == 1 ){
        return(list( c_index       = c.index(u5 = u5, cut_hazard1 = cut_hazard1),
                     selected.gene = (summary.fit[,1])[! is.na(summary.fit[,1])], 
                     coef.and.p    = as.data.frame(d) ))}
        else{return(list( c_index       = c.index(u5 = u5, cut_hazard1 = cut_hazard1),
                          selected.gene = (summary.fit[,1])[! is.na(summary.fit[,1])] ))}
      }
    }
    gc()  
  })
  
  #### Bind Result ####
  if(nfold == 1 && penalty.alpha != 0){return( list(   c_index          =  mean(unlist(cox.gene.penalty.cv.list[1,])),
                                                       selected.gene    =  cox.gene.penalty.cv.list["selected.gene",],
                                                       coef.and.p       =  cox.gene.penalty.cv.list[3,]))}
  else{return( list(   c_index          =  mean(unlist(cox.gene.penalty.cv.list[1,])),
                       selected.gene    =  cox.gene.penalty.cv.list["selected.gene",] ))}
  
}
