source('s_glmnet_lambda.R')
library(caret)
library(glmnet)
#source('s_gene_regression_fun.R')
f_lm_fit_and_predict<-function(selected_vars, target_col, input_data, in_train ){
    lm_formula=as.formula(f_p('%s~%s', target_col, paste0(selected_vars, collapse = '+')))
    lm_obj=lm(lm_formula, data = as.data.frame(input_data[in_train,]))
    pred_val=predict.lm(lm_obj, newdata = as.data.frame(input_data[-in_train,]))
    pred_val
}

f_nest_cv_glmnet <- function(input_data_raw, target_col, alpha_list = c(1), nfolds = 10, penalty_factors, my_norm = TRUE, debug = FALSE){

    if (my_norm == TRUE){
        input_data = as.matrix(input_data_raw)
    }else{
        input_data = scale(input_data_raw)
    }
        
    set.seed(889)
    #Split by balance y.
    #splits <- createFolds(input_data[,target_col], k=nfolds,returnTrain = TRUE)

    #Split by balance population
    splits = createFolds(sample_info[rownames(input_data), 'pop'], k=nfolds,returnTrain = TRUE)
    
    results <- lapply(splits, 
                      function(x, dat, target_col) {
                          holdout <- (1:nrow(dat))[-unique(x)]
                          data.frame(index = holdout, 
                                     obs = dat[holdout, target_col])
                      },
                      dat = input_data, target_col = target_col)
    #str(results)
    mods <- vector(mode = "list", length = length(splits))
    if(f_judge_debug(F)){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
        
    }

    
    
    train_perf = rep(0, times = nfolds)
    for(i in seq(along = splits)) {
        in_train <- unique(splits[[i]])
        set.seed(2)
        a=f_glmnet_max_lambda(input_data[in_train, colnames(input_data) != target_col], input_data[in_train, target_col])
        print(table(sample_info[rownames(input_data)[in_train], 'pop']))
        input_features = setdiff(colnames(input_data), c(target_col))
        test_data = input_data[-in_train,]
        train_data = input_data[in_train,]

        if (my_norm == TRUE){
            train_data_norm = scale(train_data)
            test_data_norm = scale(test_data, center = colMeans(train_data) , scale = apply(train_data, 2, sd, na.rm = TRUE))
            train_data_norm[is.na(train_data_norm)] = 0
            test_data_norm[is.na(test_data_norm)] =0
        }else{
            test_data_norm  = test_data
            train_data_norm = train_data
        }


        mod<- f_cv_glmnet(alpha_list, train_data_norm, target_col, fold_num = nfolds, penalty_factors)
        
        
        results[[i]]$pred <- as.vector(predict.cv.glmnet(mod, test_data_norm[, input_features], s='lambda.1se')[,1])
        results[[i]]$obs <- test_data_norm[,target_col]
        #predict(mod, input_data[-in_train,colnames(input_data)!=target_col ] , s='lambda.1se')
        glm_coef =coef(mod)
        
        print(f_p('fold %s lambda range(%s, %s), 1se %s', i, min(mod$lambda), max(mod$lambda), mod$lambda.1se))
        print(glm_coef[glm_coef[,1]!=0,])
        selected_vars = setdiff( names(glm_coef[glm_coef[,1]!=0,]), '(Intercept)')
        #results[[i]]$pred = f_lm_fit_and_predict(selected_vars, target_col, input_data, in_train )
        mods[[i]] <- mod
        train_perf[i] = mod$train_perf[1,'Rsquared']
    }
    #Use all the data to fit the model

    


    ##Creat a list that is in caret output format.
    fit_list = list()
    fit_list$results = f_convert_cv_results_to_SD(results, train_perf, quite = FALSE)    
    fit_list$train_perf=data.frame(Rsquared=mean(train_perf), RsquaredSD = sd(train_perf))

    scaled_input_data = scale(input_data)
    scaled_input_data[is.na(scaled_input_data)]= 0
    fit_list$finalModel = f_cv_glmnet(alpha_list, scaled_input_data, target_col, penalty_factors = penalty_factors)
    fit_list$bestTune = data.frame(alpha = fit_list$finalModel$alpha, lambda = fit_list$finalModel$lambda.1se)
    fit_list$pred = f_get_cv_prediction(fit_list$finalModel, scaled_input_data[, target_col] )
    fit_list$scale_mean = attr(scaled_input_data, 'scaled:center')
    fit_list$scale_sd = attr(scaled_input_data, 'scaled:scale')


    #Scale the feature coef back to the 
    glm_coef=as.data.frame(as.matrix(coef(fit_list$finalModel, fit_list$bestTune$lambda)))
    colnames(glm_coef) = c('coef')
    input_features = grep('Intercept', rownames(glm_coef), value = T , invert = T)
    key_features = glm_coef[,'coef']
    names(key_features) = rownames(glm_coef)
    fit_list$key_features = key_features

    return (fit_list)
}


library(stringr)
f_cv_glmnet<-function(alpha_list, train_data, target_col, fold_num = 10, penalty_factors){
    alpha_models <- vector(mode = "list", length = length(alpha_list))
    cv_performance <- rep(0, times = length(alpha_list))
    names(alpha_models) = paste0('alpha', alpha_list)
    names(cv_performance) = names(alpha_models)
    for (loc_alpha in alpha_list){
        loc_name = f_p('alpha%s', loc_alpha)
        set.seed(998)
                
        alpha_models[[loc_name]]=cv.glmnet(x = train_data[, colnames(train_data)!=target_col],
                                           y = train_data[, target_col], alpha = loc_alpha, nfolds = fold_num, nlambda = 100,
                                           keep = TRUE, penalty.factor = penalty_factors, standardize = T, maxit = 50000)#, lambda = (200:1)/100)
        cv_performance[[loc_name]] = min(alpha_models[[loc_name]]$cvm)
    }
    best_alpha = which.min(cv_performance)
    best_fit=alpha_models[[best_alpha]]
    best_fit$alpha = as.numeric(str_replace(names(cv_performance)[best_alpha], 'alpha', ''))
    
    best_fit$train_perf = f_get_glmnet_training_stats(best_fit, train_data, target_col )
    return (best_fit)    
}





   
f_get_cv_stats <- function(cv_fit, input_data, target_col){
    #Extract the CV testing performance in the cv.glmnet.
    results=list()
    for (i in 1:max(cv_fit$foldid)){
        fold_index = which(cv_fit$foldid  == i)    
        results[[f_p('fold%s', i)]] = data.frame(obs = input_data[fold_index, target_col], pred = cv_fit$fit.preval[fold_index, which(cv_fit$lambda == cv_fit$lambda.1se) ])
    }

    f_convert_cv_results_to_SD(results)
}

f_get_glmnet_training_stats <- function(cv_fit, input_data, target_col, debug = FALSE){
    #Use the cv_fit to fit the original training data. Get the training error
    if (f_judge_debug(debug)){
    ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
    }
    cv_pred <- as.vector(predict.cv.glmnet(cv_fit, input_data[,colnames(input_data)!=target_col ] , s='lambda.1se')[,1])
    result_df = data.frame(pred = cv_pred, obs = input_data[, target_col])
    f_tradition_rquare(result_df, R2_method)
} 

f_convert_cv_results_to_SD <- function(results, train_perf,quite = TRUE){
    library(plyr)
    library(dplyr)
        
    stats_df = ldply(lapply(results, f_tradition_rquare, R2_method ))
    stats_df$rsquared = stats_df$Rsquared
    stats_df$train = train_perf
    stats_df=stats_df[,c('Rsquared', 'train', 'RSME', 'rsquared')]
    if (quite == FALSE){
        print (stats_df)
    }
    stats_df %>% dplyr::summarise(Rsquared = mean(rsquared), RsquaredSD = sd(rsquared), train = mean(train), trainSD = sd(train) )
}

f_get_cv_prediction <-function(cv_fit, obs){
    pred = cv_fit$fit.preval[, which(cv_fit$lambda == cv_fit$lambda.1se)]
    pred_df = data.frame(pred = pred, obs = obs, rowIndex = 1:length(obs), fold = cv_fit$foldid)
    return (pred_df)
}

f_tradition_rquare<- function(data, method = 'traditional'){
    #Move to corr because the validation data are in different scales of the train data.
    if (sd(data$pred) == 0 & method == 'corr' ){
        r_square = 0
    }else{
        r_square = R2(data$pred, data$obs, formula = method, na.rm = FALSE)
    }
    data.frame( Rsquared=r_square, RSME = RMSE(data$pred, data$obs, na.rm = FALSE))
}

f_tradition_rquare2<- function(pred, obs, method = 'traditional'){
    #Move to corr because the validation data are in different scales of the train data.
        
    if (sd(pred) == 0 & method == 'corr'){
        r_square = 0
    }else{
        r_square = R2(pred, obs, formula = method, na.rm = FALSE)
    }
    r_square
}



