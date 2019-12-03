#######################################################################
# Logistic Model 
#######################################################################
library(glmnet)
library(data.table)
library(selectiveInference)

#################################################

lr.gene.penalty.cv.result = function(user.data, response.col, gene.type, 
                                     padj.method   = "fdr", fdr.threshold  = 0.05,
                                     penalty.alpha = 1 , nfold = 1, force.add.cov = NULL ){
  set.seed(2)
  t       = sample(1:nfold,size = nrow(user.data), replace = T)
  gene.col             = gene.type$gene.var
  user.data[,gene.col] = apply(user.data[,gene.col], 2, as.numeric)
  response             = as.factor(user.data[,response.col]) 
  pred.result          = rep(0, nrow(user.data) )
  adaptive.lasso       = (penalty.alpha == 2)
  gc()
  
  #### Training/Testing ####
  lr.gene.penalty.cv.list = sapply(1:nfold, function(fold){
    train  = !(t==fold)
    test   =  (t==fold)
    if(nfold==1){train  = (t==fold);test   =  (t==fold)} # fold = 1
    
    #### Delete the 0 sum gene covariate column ####
    tmpsum = apply(user.data[ train, gene.col ], 2, sum)
    gene.col.train = gene.col[tmpsum != 0]
    
    #### univariate logistic regression fdr ####
    if( padj.method == "fdr" ){
      tmp.p = rep(NA, length(gene.col.train))
      tmp.p = sapply( gene.col.train, function(i){ summary( glm( response ~ . , family=binomial, data = data.frame(response = as.numeric(response == 1), as.numeric(user.data[,i])), subset = train ) )$coefficients[2,4] })
      tmp.p = p.adjust(tmp.p, padj.method) 
      fdr.genes.col = gene.col.train[ which( tmp.p <= fdr.threshold )]
      if (length(fdr.genes.col) == 0) stop("small fdr.threshold --> no genes selected")}
    else{
        fdr.genes.col = gene.col.train
      }
    
    #### Regression part ####
    x <- model.matrix( response ~ . , data.frame( response = as.numeric(response == 1), user.data[,fdr.genes.col]) )[,-1]
    y <- as.numeric(response == 1)
    
    x.train = x[train,]
    x.test  = x[test ,]
    
    y.train = y[train]
    y.test  = y[test ]
    
    if( adaptive.lasso ){penalty.alpha = 1}
    
    set.seed(6)
    cv.out       <- cv.glmnet(x.train, y.train, alpha = penalty.alpha, family="binomial", type.measure = "mse", nfolds = 5 )
    cv.est.coef       = coef(cv.out, s = cv.out$lambda.min)
    best_coef         = as.numeric(cv.est.coef)
    geneSel           = rownames(as.matrix(cv.est.coef))[which(!( best_coef == 0))][-1]
    gc()
    lambda_min   <- cv.out$lambda.min
    if(sum(as.numeric(coef(cv.out, s = lambda_min))==0) == length(fdr.genes.col)){
      lambda_min <- cv.out$lambda[as.numeric(which(cv.out$nzero != 0)[1])]
      }
    
    #### Ridge and Lasso ####
    if( ! adaptive.lasso ){
      
      #### Standard error and p-value for estimated coefficients ####
      if( nfold == 1 && penalty.alpha == 1 ){
        cv.est.coef       = coef(cv.out, s = lambda_min)
        best_coef         = as.numeric(cv.est.coef)
        #li = fixedLassoInf(x = x, y = y, beta = best_coef, lambda = cv.out$lambda.min, family = "binomial")
        #li = fixedLassoInf(x = x[,-which(duplicated(x, MARGIN = 2))], y = y, 
        #                  beta = best_coef[-which(duplicated(x, MARGIN = 2))], 
        #                  lambda = lambda_min, family = "binomial")
        best_coef1 = best_coef[best_coef!=0]
        li = fixedLassoInf(x = x[,colnames(x) %in% geneSel], y = y,
                           beta = best_coef1,#cv.est.coef,#best_coef, 
                           lambda = lambda_min, family = "binomial")
        best_coef_p_value = li$pv
        
        la_geneSel        = row.names(cv.est.coef)[which(best_coef != 0)][-1]
        d = cbind(gene_list = la_geneSel, 
                  estimated_coefficient = best_coef[which(best_coef != 0)][-1], 
                  p_value = best_coef_p_value)
        }
      
      #### training data: find cut point for testing data ####
      lasso_prob = predict(cv.out , newx = x.train , s = lambda_min , type = "response")
      
      #### 0.01-0.99 threshold find cut point by Youden index ####
      training_sen = rep(0, 9999)
      training_spe = rep(0, 9999)
      for(i in 1:9999){
        training_response_by_threshold = rep(0, length(lasso_prob))
        over_threshold_index = which(lasso_prob > 0.0001*i)
        training_response_by_threshold[over_threshold_index] = 1
        
        ## Confusion Table
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

      
      #### testing data: use cut point find (1)sen spe accuracy and (2)testing AUC ####
      pred.prob         = predict(cv.out , newx = x.test , s=cv.out$lambda.min, type="response") 
      
      #### use training cut point to find sen spe accuracy ####
      testing_response_by_cutpt = rep(0, length(pred.prob))
      over_cutpt_index = which(pred.prob > cutpt)
      testing_response_by_cutpt[over_cutpt_index] = 1
      
      testing_sen_result = sum((testing_response_by_cutpt == 1)*(y.test == 1))/sum(y.test == 1)
      testing_spe_result = sum((testing_response_by_cutpt == 0)*(y.test == 0))/sum(y.test == 0)
      
      #### 0.01-0.99 threshold find ROC curve ####
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
      
      if(is.null(force.add.cov)){
        if( nfold == 1 && penalty.alpha == 1 ){
        return(list( cp = cutpt,
                     testing_pre0_res0 = sum((testing_response_by_cutpt == 0)*(y.test == 0)),
                     testing_pre0_res1 = sum((testing_response_by_cutpt == 0)*(y.test == 1)),
                     testing_pre1_res0 = sum((testing_response_by_cutpt == 1)*(y.test == 0)),
                     testing_pre1_res1 = sum((testing_response_by_cutpt == 1)*(y.test == 1)),
                     test.sen      = testing_sen_result,
                     test.spe      = testing_spe_result,
                     test.auc      = as.numeric(rocobj$auc),
                     selected.gene = ( (as.matrix(coef(cv.out, s = lambda_min)))[-1,] )[which( (as.numeric(coef(cv.out, s = lambda_min))[-1]) != 0)],
                     coef.and.p    = as.data.frame(d) 
                     ) )
        }else{
          return(list( cp = cutpt,
                       testing_pre0_res0 = sum((testing_response_by_cutpt == 0)*(y.test == 0)),
                       testing_pre0_res1 = sum((testing_response_by_cutpt == 0)*(y.test == 1)),
                       testing_pre1_res0 = sum((testing_response_by_cutpt == 1)*(y.test == 0)),
                       testing_pre1_res1 = sum((testing_response_by_cutpt == 1)*(y.test == 1)),
                       test.sen      = testing_sen_result,
                       test.spe      = testing_spe_result,
                       test.auc      = as.numeric(rocobj$auc),
                       selected.gene = ( (as.matrix(coef(cv.out, s = lambda_min)))[-1,] )[which( (as.numeric(coef(cv.out, s = lambda_min))[-1]) != 0)]
          ) )
        }
      }else{
        #### training data: find cut point for testing data ####
        geneSel              = names(( (as.matrix(coef(cv.out, s = lambda_min)))[-1,] )[which( (as.numeric(coef(cv.out, s = lambda_min))[-1]) != 0)])
        lr.data              = data.frame(cbind(y.train = y.train, user.data[train, force.add.cov], x.train[, c(which(colnames(x.train) %in% geneSel)) ]) )
        fit.geneSel.forceCov = glm( y.train ~ . , family = "binomial", data = lr.data,  control = list(maxit = 20000))
        lr_prob              = predict(fit.geneSel.forceCov, newdata =  data.frame(cbind(y.train = y.train, user.data[train, force.add.cov], x.train[, c(which(colnames(x.train) %in% geneSel)) ]) ), type="response")
        
        #### 0.01-0.99 threshold find cut point by Youden index ####
        training_sen = rep(0, 9999)
        training_spe = rep(0, 9999)
        for(i in 1:9999){
          training_response_by_threshold = rep(0, length(lr_prob))
          over_threshold_index = which(lasso_prob > 0.0001*i)
          training_response_by_threshold[over_threshold_index] = 1
          
          ## Confusion Table
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
        
        #### testing data: use cut point find (1)sen spe accuracy and (2)testing AUC ####
        lr_prob              = predict(fit.geneSel.forceCov, newdata =  data.frame(cbind(y.train = y.test , user.data[test , force.add.cov],x.test [, c(which(colnames(x.test ) %in% geneSel)) ]) ), type="response")
        
        #### use training cut point to find sen spe accuracy ####
        testing_response_by_cutpt = rep(0, length(lr_prob))
        over_cutpt_index = which(lr_prob > cutpt)
        testing_response_by_cutpt[over_cutpt_index] = 1
        
        testing_sen_result = sum((testing_response_by_cutpt == 1)*(y.test == 1))/sum(y.test == 1)
        testing_spe_result = sum((testing_response_by_cutpt == 0)*(y.test == 0))/sum(y.test == 0)
        
        #### 0.01-0.99 threshold find ROC curve ####
        testing_sen = rep(0, 9999)
        testing_spe = rep(0, 9999)
        for(i in 1:9999){
          testing_response_by_threshold = rep(0, length(lr_prob))
          over_threshold_index = which(lr_prob > 0.0001*i)
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
        if( nfold == 1 && penalty.alpha == 1 ){
        return(list( cp = cutpt,
                     testing_pre0_res0 = sum((testing_response_by_cutpt == 0)*(y.test == 0)),
                     testing_pre0_res1 = sum((testing_response_by_cutpt == 0)*(y.test == 1)),
                     testing_pre1_res0 = sum((testing_response_by_cutpt == 1)*(y.test == 0)),
                     testing_pre1_res1 = sum((testing_response_by_cutpt == 1)*(y.test == 1)),
                     test.sen      = testing_sen_result,
                     test.spe      = testing_spe_result,
                     test.auc      = as.numeric(rocobj$auc),
                     selected.gene = summary(fit.geneSel.forceCov)$coefficients[,1],
                     coef.and.p    = as.data.frame(d) ))
        }else{
          return(list( cp = cutpt,
                       testing_pre0_res0 = sum((testing_response_by_cutpt == 0)*(y.test == 0)),
                       testing_pre0_res1 = sum((testing_response_by_cutpt == 0)*(y.test == 1)),
                       testing_pre1_res0 = sum((testing_response_by_cutpt == 1)*(y.test == 0)),
                       testing_pre1_res1 = sum((testing_response_by_cutpt == 1)*(y.test == 1)),
                       test.sen      = testing_sen_result,
                       test.spe      = testing_spe_result,
                       test.auc      = as.numeric(rocobj$auc),
                       selected.gene = summary(fit.geneSel.forceCov)$coefficients[,1]
                       ))
        }
      }
      
    }else{
      #### Adaptive Lasso ####
      #### Estimated coefficients ####
      cv.est.coef       = coef(cv.out, s = cv.out$lambda.min)
      best_coef         = as.numeric(cv.est.coef)
      geneSel           = rownames(as.matrix(cv.est.coef))[which(!( best_coef == 0))][-1]
      
      #### Weight Change ####
      w.vector          = 1 / (abs(best_coef)[which(best_coef != 0)][-1])
      w.vector[w.vector == Inf] = 999999999
      
      #### Perform adaptive LASSO ####
      fit1_cv                   = cv.glmnet(x = x.train[,colnames(x.train) %in% geneSel], y = y.train, family = "binomial", nfold = 5, alpha = 1, penalty.factor = w.vector)
      lambda_min                = fit1_cv$lambda.min
      cv.adap.lasso.est.coef    = coef(fit1_cv, s = fit1_cv$lambda.min)
      best_coef_adap_la         = as.numeric(cv.adap.lasso.est.coef)
      adap_la_geneSel           = rownames(as.matrix(cv.adap.lasso.est.coef))[which(!( best_coef_adap_la == 0))][-1]
      
      #### Standard error and p-value for estimated coefficients ####
      if( nfold == 1 && penalty.alpha == 1 ){
        adap.lasso.est.coef       = as.numeric(coef(fit1_cv, s = fit1_cv$lambda.min))
        
        cv.adap.lasso.est.coef    = coef(fit1_cv, s = fit1_cv$lambda.min)
        best_coef_adap_la         = as.numeric(cv.adap.lasso.est.coef)
        ###
        #cv.adap.lasso.est.coef1 = cv.adap.lasso.est.coef[cv.adap.lasso.est.coef!=0]
        ###
        adap.l.i                  = fixedLassoInf(x = x[,colnames(x) %in% geneSel], y = y, beta = cv.adap.lasso.est.coef, lambda = fit1_cv$lambda.min, family = "binomial")
        best_coef_adap_la_p_value = adap.l.i$pv
        
        la_geneSel          = row.names(cv.adap.lasso.est.coef)[which(best_coef_adap_la != 0)][-1]
        d = cbind(gene_list = la_geneSel,
                  estimated_coefficient = best_coef_adap_la[which(best_coef_adap_la != 0)][-1], 
                  p_value = best_coef_adap_la_p_value)
        }
      
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
      if(is.null(force.add.cov)){
        if( nfold == 1 && penalty.alpha == 1 ){
        return(list( cp = cutpt,
                     testing_pre0_res0 = sum((testing_response_by_cutpt == 0)*(y.test == 0)),
                     testing_pre0_res1 = sum((testing_response_by_cutpt == 0)*(y.test == 1)),
                     testing_pre1_res0 = sum((testing_response_by_cutpt == 1)*(y.test == 0)),
                     testing_pre1_res1 = sum((testing_response_by_cutpt == 1)*(y.test == 1)),
                     test.sen      = testing_sen_result,
                     test.spe      = testing_spe_result,
                     test.auc      = as.numeric(rocobj$auc),
                     selected.gene = (as.matrix(cv.adap.lasso.est.coef))[which(cv.adap.lasso.est.coef != 0),],
                     coef.and.p    = as.data.frame(d)
                     ))
        }else{
          return(list( cp = cutpt,
                       testing_pre0_res0 = sum((testing_response_by_cutpt == 0)*(y.test == 0)),
                       testing_pre0_res1 = sum((testing_response_by_cutpt == 0)*(y.test == 1)),
                       testing_pre1_res0 = sum((testing_response_by_cutpt == 1)*(y.test == 0)),
                       testing_pre1_res1 = sum((testing_response_by_cutpt == 1)*(y.test == 1)),
                       test.sen      = testing_sen_result,
                       test.spe      = testing_spe_result,
                       test.auc      = as.numeric(rocobj$auc),
                       selected.gene = ( (as.matrix(coef(fit1_cv, s = lambda_min)))[-1,] )[which((as.numeric(coef(fit1_cv, s = lambda_min))[-1]) != 0) ]
          ))
        }
        
      }else{
        geneSel              = names(( (as.matrix(coef(fit1_cv, s = lambda_min)))[-1,] )[which( (as.numeric(coef(fit1_cv, s = lambda_min))[-1]) != 0)])
        lr.data              = data.frame(cbind(y.train = y.train, user.data[train, force.add.cov],x.train[, c(which(colnames(x.train) %in% geneSel)) ]) )
        fit.geneSel.forceCov = glm( y.train ~ . , family = "binomial", data = lr.data,  control = list(maxit = 20000))
        lr_prob              = predict(fit.geneSel.forceCov, newdata =  data.frame(cbind(y.train = y.train, user.data[train, force.add.cov],x.train[, c(which(colnames(x.train) %in% geneSel)) ]) ), type="response")

        #### 0.01-0.99 threshold find cut point by Youden index ####
        training_sen = rep(0, 9999)
        training_spe = rep(0, 9999)
        for(i in 1:9999){
          training_response_by_threshold = rep(0, length(lr_prob))
          over_threshold_index = which(lasso_prob > 0.0001*i)
          training_response_by_threshold[over_threshold_index] = 1
          
          ## Confusion Table
          training_sen_index = sum((training_response_by_threshold == 1)*(y.train == 1))
          training_spe_index = sum((training_response_by_threshold == 0)*(y.train == 0))
          training_sen[i] = training_sen_index/sum(y.train == 1)
          training_spe[i] = training_spe_index/sum(y.train == 0)
        }
        ## Youden index
        Youden_index = rep(0, 9999)
        for(i in 1:9999){
          Youden_index[i] = training_sen[i] + training_spe[i] - 1
        }
        ## Cut Point By Youden index
        cutpt_index = which(Youden_index == max(Youden_index))
        cutpt       = 0.0001*median(cutpt_index) # max(cutpt_index)
        
        #### testing data ####
        lr_prob              = predict(fit.geneSel.forceCov, newdata =  data.frame(cbind(y.train = y.test , user.data[test , force.add.cov],x.test [, c(which(colnames(x.test ) %in% geneSel)) ]) ), type="response")
        
        #### use training cut point to find sen spe accuracy ####
        testing_response_by_cutpt = rep(0, length(lr_prob))
        over_cutpt_index = which(lr_prob > cutpt)
        testing_response_by_cutpt[over_cutpt_index] = 1
        
        testing_sen_result = sum((testing_response_by_cutpt == 1)*(y.test == 1))/sum(y.test == 1)
        testing_spe_result = sum((testing_response_by_cutpt == 0)*(y.test == 0))/sum(y.test == 0)
        
        #### 0.01-0.99 threshold find ROC curve ####
        testing_sen = rep(0, 9999)
        testing_spe = rep(0, 9999)
        for(i in 1:9999){
          testing_response_by_threshold = rep(0, length(lr_prob))
          over_threshold_index = which(lr_prob > 0.0001*i)
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
        if( nfold == 1 && penalty.alpha == 1 ){
        return(list( cp = cutpt,
                     testing_pre0_res0 = sum((testing_response_by_cutpt == 0)*(y.test == 0)),
                     testing_pre0_res1 = sum((testing_response_by_cutpt == 0)*(y.test == 1)),
                     testing_pre1_res0 = sum((testing_response_by_cutpt == 1)*(y.test == 0)),
                     testing_pre1_res1 = sum((testing_response_by_cutpt == 1)*(y.test == 1)),
                     test.sen      = testing_sen_result,
                     test.spe      = testing_spe_result,
                     test.auc      = as.numeric(rocobj$auc),
                     selected.gene = summary(fit.geneSel.forceCov)$coefficients[,1],
                     coef.and.p    = as.data.frame(d)
                     ))
        }else{
          return(list( cp = cutpt,
                       testing_pre0_res0 = sum((testing_response_by_cutpt == 0)*(y.test == 0)),
                       testing_pre0_res1 = sum((testing_response_by_cutpt == 0)*(y.test == 1)),
                       testing_pre1_res0 = sum((testing_response_by_cutpt == 1)*(y.test == 0)),
                       testing_pre1_res1 = sum((testing_response_by_cutpt == 1)*(y.test == 1)),
                       test.sen      = testing_sen_result,
                       test.spe      = testing_spe_result,
                       test.auc      = as.numeric(rocobj$auc),
                       selected.gene = summary(fit.geneSel.forceCov)$coefficients[,1]
          ))
        }
      }
    }
    
    gc()
    })
  
  ### combine the final result
  # for(i in 1:nfold){ pred.result = pred.result + unlist(lr.gene.penalty.cv.list["pred.result"])}
  cp    = unlist(lr.gene.penalty.cv.list["cp",])
  
  t_1_1 = sum(unlist(lr.gene.penalty.cv.list[2,]))
  t_1_2 = sum(unlist(lr.gene.penalty.cv.list[4,]))
  t_2_1 = sum(unlist(lr.gene.penalty.cv.list[3,]))
  t_2_2 = sum(unlist(lr.gene.penalty.cv.list[5,]))
  tmp = matrix(c(t_1_1,t_2_1,t_1_2,t_2_2),2,2)
  
  if(nfold == 1 && penalty.alpha != 0){return( list( cutpoint      =  cp,
                                                     sen           =  tmp[2,2]/sum(tmp[2,]),
                                                     spe           =  tmp[1,1]/sum(tmp[1,]),
                                                     accuracy      =  sum(diag(tmp))/sum(tmp),
                                                     mean.auc      =  mean(unlist(lr.gene.penalty.cv.list[8,])),
                                                     selected.gene =  lr.gene.penalty.cv.list[9,],
                                                     coef.and.p    =  lr.gene.penalty.cv.list[10,])
                                               )}
  else{
    return( list( cutpoint      =  cp,
                  sen           =  tmp[2,2]/sum(tmp[2,]),
                  spe           =  tmp[1,1]/sum(tmp[1,]),
                  accuracy      =  sum(diag(tmp))/sum(tmp),
                  mean.auc      =  mean(unlist(lr.gene.penalty.cv.list[8,])),
                  selected.gene =  lr.gene.penalty.cv.list[9,]
    ))
  }
  
}
