library(randomForest)
library(mclust)

#my_train = final_train_data
#target_col = 'gene.RNASEQ'
f_get_server_name <- function(){
    return (Sys.info()["nodename"])
}



f_caret_regression_return_fit <- function(my_train, target_col, learner_name, tuneGrid,quite=FALSE)
{
  
  #cl <- makeCluster(3)
  #registerDoSNOW(cl)
#    registerDoMC(cores = 2)
    library(caret)
    cat('learner name:', learner_name, '\n')
    dim(my_train)
    head(my_train[,1:10])
    zero_cols = nearZeroVar(my_train, freqCut = 90/10)
    if (length(zero_cols) == 0)
        non_zero_data = my_train
    else
        non_zero_data = my_train[, -zero_cols]
    
    rm(my_train)
    
    #numeric_cols = colnames(non_zero_data)[sapply(non_zero_data, is.numeric)]
    #correlations <- cor(non_zero_data[, numeric_cols])
    #head(correlations[,1:10])
    #other_cols = setdiff(colnames(non_zero_data), numeric_cols)
    #highCorr <- findCorrelation(correlations, cutoff = .9)

    
    
    #if(length(highCorr) > 0){
    #    cat('Correlated columns:', length(highCorr), 'out of', ncol(non_zero_data), '\n')
    #    remain_cols = setdiff(colnames(non_zero_data), colnames(correlations)[highCorr])
    #    #non_zero_data = non_zero_data[,  unique(c(remain_cols, target_col))]
    #}
    
    #non_zero_data = my_train
    
    #dummies <- dummyVars(as.formula(f_p('%s ~ .', target_col)), data = non_zero_data1)
    #non_zero_data=as.data.frame(predict(dummies, newdata = non_zero_data1))

    
    cat('After filtering, training size', dim(non_zero_data), '(features, samples)')
    #colnames(non_zero_data)

    #if (f_get_server_name() == 'loire'){
    #    registerDoMC(cores = 2)
    #}else{
    #    core_num = 1 #min(6, round(20000/ncol(non_zero_data)))
    #    #registerDoMC(cores = core_num - 2)
    #    cat('Core number', core_num, '\n')
    #}
        
   
 
  fitControl <- trainControl(method = "repeatedcv",
                             number = 5,
                             repeats = 3,
                             allowParallel= TRUE)
  
  if (learner_name == 'enet')
  {
    
    Fit <- train(form = as.formula(f_p('%s ~ .', target_col)), data = non_zero_data,
                 method = "enet",
                 tuneGrid = tuneGrid,
                 metric = 'Rsquared',
                 maximize = TRUE,
                 trControl = fitControl,
                 preProc = c("center", "scale"))
    if (quite == FALSE){
        print(Fit$resample)
        #summary(Fit)
    }


    ###Get the coefficient from the model########
    
    coeff=predict(Fit$finalModel, type='coefficient', s=Fit$bestTune$fraction, mode='f')
    key_features = (coeff$coefficients[coeff$coefficients != 0])

    #reguired output for the fit.
    Fit$key_features = key_features
    Fit$max_performance=Fit$results[rownames(Fit$bestTune), 'Rsquared']
    Fit$num_features = ncol(non_zero_data)
    Fit$selected_features = length(key_features)
    #train_pred = predict(Fit, newdata = my_train, type='raw')
    
    #processed=preProcess(my_train, method = c('scale','center'))
    #processed_data=predict(processed, my_train)
    #result= data.frame(truth  = my_train[[target_col]])
    #result[,predictors(Fit)] =   processed_data[,predictors(Fit)]
    
  }else if(learner_name == 'rfe'){
      
      ctrl <- rfeControl(method = "cv", number = 5,
                       verbose = FALSE,
                       functions = rfFuncs)
     
      predVars = setdiff(colnames(non_zero_data),c(target_col))
      varSeq <- seq(2, length(predVars)-1, by = 10)
      rfRFE <- rfe(x = non_zero_data[, predVars], y = non_zero_data[,target_col], sizes = varSeq, metric = "Rsquared", rfeControl = ctrl,
                   ## now pass options to randomForest()
                   ntree = 2000)

      Fit = list()
      Fit$key_features=(rfRFE$fit$importance[,'IncNodePurity'])
      Fit$max_performance = max(rfRFE$results$Rsquared)
      
  }else if(learner_name =='glmnet'){
    Fit <- caret::train( as.formula(f_p('%s ~ .', target_col)), data = non_zero_data,
                        metric = "RMSE", 
                         method = learner_name,
                         trControl = fitControl,
                        preProc = c("center", "scale"))

    
    Fit$max_performance=Fit$results[rownames(Fit$bestTune), 'Rsquared']
    Fit$num_features = ncol(non_zero_data)
    
    glm_coef=as.data.frame(as.matrix(coef(Fit$finalModel, Fit$bestTune$lambda)))
    colnames(glm_coef) = c('coef')
    sig_coef=subset(glm_coef, coef!=0)
    key_features = sig_coef[,'coef']
    names(key_features) = rownames(sig_coef)
    Fit$key_features = key_features

    Fit$num_features = ncol(non_zero_data)
    Fit$selected_features = length(key_features)
    #rfFit = randomForest(as.formula(f_p('%s ~ .', target_col)), data = my_train, ntree = 1000, do.trace=100)
    
    #train_pred = rfFit$pred
    #result= data.frame(truth  = my_train[[target_col]])
  }else{
      cat('unknown learner name:', learner_name, '\n')
  }
  
  
  
  return (Fit)
  
}



#target_col='mixture_class'
f_caret_classification_return_fit <- function(my_train, target_col, learner_name, tuneGrid,quite=FALSE, metric = 'Accuracy')
{
    library(caret)
    library(glmnet)
    cat('learner name:', learner_name, '\n')
    dim(my_train)
    head(my_train[,1:10])
    zero_cols = nearZeroVar(my_train, freqCut = 90/10)
    if (length(zero_cols) == 0)
        non_zero_data = my_train
    else
        non_zero_data = my_train[, -zero_cols]
    
    #rm(my_train)
    
    numeric_cols = colnames(non_zero_data)[sapply(non_zero_data, is.numeric)]
    correlations <- cor(non_zero_data[, numeric_cols])
    #head(correlations[,1:10])
    other_cols = setdiff(colnames(non_zero_data), numeric_cols)
    highCorr <- findCorrelation(correlations, cutoff = .9)

    
    
    if(length(highCorr) > 0){
        cat('Correlated columns:', length(highCorr), 'out of', ncol(non_zero_data), '\n')
        remain_cols = setdiff(colnames(non_zero_data), colnames(correlations)[highCorr])
        #non_zero_data = non_zero_data[,  unique(c(remain_cols, target_col))]
    }
    
    #non_zero_data = my_train
    
    #dummies <- dummyVars(as.formula(f_p('%s ~ .', target_col)), data = non_zero_data1)
    #non_zero_data=as.data.frame(predict(dummies, newdata = non_zero_data1))

    
    cat('After filtering, training size', dim(non_zero_data), '(features, samples)')
    #colnames(non_zero_data)

    #if (f_get_server_name() == 'loire'){
    #    registerDoMC(cores = 2)
    #}else{
    #    core_num = 1 #min(6, round(20000/ncol(non_zero_data)))
    #    #registerDoMC(cores = core_num - 2)
    #    cat('Core number', core_num, '\n')
    #}
        
   
  #library(caret)
  fitControl <- trainControl(method = "repeatedcv",
                             number = 5,
                             repeats = 1,
                             allowParallel= TRUE)
  
   
  if(learner_name == 'rfe'){
      
      ctrl <- rfeControl(method = "cv", number = 5,
                       verbose = TRUE,
                       functions = rfFuncs)
     
      predVars = setdiff(colnames(non_zero_data),c(target_col))
      varSeq <- seq(2, length(predVars)-1, by = 20)
      rfRFE <- rfe(x = non_zero_data[, predVars], y = non_zero_data[,target_col], sizes = varSeq, metric = 'Accuracy', rfeControl = ctrl,
                   ## now pass options to randomForest()
                   ntree = 1000)

      rfRFE$results
      f_sort_by_col(rfRFE$fit$importance, index ='MeanDecreaseGini' )
      Fit = list()
      Fit$fit = rfRFE
      Fit$key_features=(rfRFE$fit$importance[,'MeanDecreaseGini'])
      Fit$max_performance = max(rfRFE$results$Accuracy)
      
  }else{
    Fit <- caret::train( as.formula(f_p('%s ~ .', target_col)), data = non_zero_data,
                        metric = "RMSE", 
                         method = learner_name,
                         trControl = fitControl,
                         verbose = TRUE,
                        preProc = c("center", "scale"))
    
    #rfFit = randomForest(as.formula(f_p('%s ~ .', target_col)), data = my_train, ntree = 1000, do.trace=100)
    
    #train_pred = rfFit$pred
    #result= data.frame(truth  = my_train[[target_col]])
  }
  
  
  
  return (Fit)
  
}










f_get_ecdf_value <- function(loc_array){
    return (ecdf(loc_array)(loc_array))
}

f_one_pair_tf_interaction <- function(match_line, sample_cols, tf_regions){

    #print(match_line)
    line1 = unlist(tf_regions[match_line[1],sample_cols])
    line2 = unlist(tf_regions[match_line[2],sample_cols])

    new_line= tf_regions[match_line[1],]
    new_line[, sample_cols] = f_get_ecdf_value(line1) * f_get_ecdf_value(line2)
    new_line[,'feature'] = paste0(tf_regions[match_line[1],'feature_tf'], '-' ,tf_regions[match_line[2],'feature_tf'])
    new_line[,'type'] = 'TF-TF'
    return (new_line)
}


f_convet_opts_to_output_dir <- function(opt){
    
    useful_opts = grep('(batch_name|help|test|gene|chr|output_mode)',names(opt), value = TRUE, invert = TRUE)
    useful_df=ldply(opt[useful_opts])
    useful_df$rename=str_replace_all(useful_df$.id, pattern = '_', replacement = '.')
    useful_df[useful_df$V1 == 'FALSE','rename'] = str_replace_all(useful_df[useful_df$V1 == 'FALSE','rename'], 'add', 'rm')
    other_opts = grep('TRUE|FALSE', useful_df$V1, invert = TRUE)
    useful_df[other_opts, 'rename'] = paste0(useful_df[other_opts, 'rename'], '.', useful_df[other_opts, 'V1'])
    #str(useful_df)
    dir_name = paste(useful_df$rename, collapse = '_')
    return (dir_name)
}



