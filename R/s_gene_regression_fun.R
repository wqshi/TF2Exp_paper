library(randomForest, quietly = TRUE)
library(mclust, quietly = TRUE)
library(caret, quietly = TRUE)
source('~/R/s_function.R', chdir = TRUE)
library(GenomicRanges)

tuneGridList = list(
     gbm = expand.grid(.interaction.depth = seq(1, 7, by = 2),
                            .n.trees = seq(100, 1000, by = 50),
                            .shrinkage = c(0.01, 0.1)),
    enet = expand.grid(.lambda = c(0, 0.01, .1),
                            .fraction = seq(.05, 1, length = 20)),

    glmnet = expand.grid(.alpha = c(0, .1, .2, .4, .6, .8, 1),
                          .lambda = c(0.0001, 0.001, 0.005 ,seq(.01, .2, length = 20)))
)

f_get_server_name <- function(){
    return (Sys.info()["nodename"])
}


#my_train = final_train_data
#target_col = 'gene.RNASEQ'
#learner_name = 'rf'
f_caret_regression_return_fit <- function(my_train, target_col, learner_name, tuneGrid,quite=FALSE)
{
    
  #cl <- makeCluster(3)
  #registerDoSNOW(cl)
#    registerDoMC(cores = 2)
    library(caret, quietly = TRUE)
    #set.seed(3456)
    cat('learner name:', learner_name, '\n')
    dim(my_train)
    head(my_train[,1:10])
    zero_cols = nearZeroVar(my_train, freqCut = 90/10)
    if (length(zero_cols) == 0)
        non_zero_data = my_train
    else
        non_zero_data = my_train[, -zero_cols]

    if (!target_col %in% colnames(non_zero_data)){
        non_zero_data[target_col] = my_train[target_col]
    }
    
    rm(my_train)
    
    rm_high_cor = FALSE
    if (rm_high_cor == TRUE){
        numeric_cols = colnames(non_zero_data)[sapply(non_zero_data, is.numeric)]
        correlations <- cor(non_zero_data[, numeric_cols])
        head(correlations[,1:10])
        other_cols = setdiff(colnames(non_zero_data), numeric_cols)
        highCorr <- findCorrelation(correlations, cutoff = .9)

    
    
        if(length(highCorr) > 0){
            cat('Correlated columns:', length(highCorr), 'out of', ncol(non_zero_data), '\n')
            remain_cols = setdiff(colnames(non_zero_data), colnames(correlations)[highCorr])
            #non_zero_data = non_zero_data[,  unique(c(remain_cols, target_col))]
        }
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
        
   
 
  fitControl <- trainControl(method = "repeatedcv",
                             number = 5,
                             repeats = 1,
                             savePredictions = 'final',
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
    
  }else if(learner_name == 'rf'){    
      ctrl <- rfeControl(method = "cv", number = 5,
                       verbose = FALSE,
                       functions = rfFuncs)
     
      predVars = setdiff(colnames(non_zero_data),c(target_col))
      varSeq <- seq(2, length(predVars)-1, length.out = 5)
      rfRFE <- rfe(x = non_zero_data[, predVars], y = non_zero_data[,target_col],
                   sizes = varSeq, metric = "Rsquared", rfeControl = ctrl,
                   ## now pass options to randomForest()
                   ntree = 1000,
                   trControl = fitControl,
                   tuneLength = 6)

      Fit = rfRFE
      Fit$bestTune=Fit$results[Fit$results$Variables == Fit$bestSubset,]
      
      Fit$pred=data.frame(rowIndex = 1:nrow(non_zero_data), pred = rfRFE$fit$predicted, obs = non_zero_data[, target_col]  )          
      Fit$key_features=(rfRFE$fit$importance[,'IncNodePurity'])
      Fit$max_performance = max(rfRFE$results$Rsquared)
  }else if(learner_name == 'svm'){
      #Not work
      ctrl <- rfeControl(method = "cv", number = 5,
                       verbose = FALSE,
                       functions = caretFuncs)
     
      predVars = setdiff(colnames(non_zero_data),c(target_col))
      varSeq <- seq(2, length(predVars)-1, by = 20)
      rfRFE <- rfe(x = non_zero_data[, predVars], y = non_zero_data[,target_col],
                   sizes = varSeq, metric = "Rsquared", rfeControl = ctrl,
                   method = "svmRadial",
                   tuneLength = 12
                   )

      Fit = list()
      Fit$key_features=(rfRFE$fit$importance[,'IncNodePurity'])
      Fit$max_performance = max(rfRFE$results$Rsquared)
  }else if(learner_name =='glmnet'){
      
    Fit <- caret::train( as.formula(f_p('%s ~ .', target_col)), data = non_zero_data,
                        metric = "RMSE", 
                         method = 'glmnet',
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


  }else if(learner_name == 'gbm'){
     
    Fit <- caret::train( as.formula(f_p('%s ~ .', target_col)), data = non_zero_data,
                        metric = "RMSE", 
                        method = 'gbm',
                        tuneGrid = tuneGrid,
                        trControl = fitControl,
                        preProc = c("center", "scale"))

    Fit$results
    gbm_coef = summary(Fit$finalModel, plotit = FALSE )
    
    Fit$max_performance=Fit$results[rownames(Fit$bestTune), 'Rsquared']
    Fit$num_features = ncol(non_zero_data)

    head(gbm_coef)
    colnames(gbm_coef) = c('var','coef')
    sig_coef=subset(gbm_coef, coef!=0)
    key_features = sig_coef[,'coef']
    names(key_features) = rownames(sig_coef)
    Fit$key_features = key_features  
           
  }else{
      cat('unknown learner name:', learner_name, '\n')
  }
  
  Fit$num_features = ncol(non_zero_data)
  Fit$selected_features = length(Fit$key_features)
  
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
                             number = 10,
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
    
    useful_opts = grep('(batch_name|help|test|gene|chr|output_mode|permutation|add_TF_exp|add_predict_tf)',names(opt), value = TRUE, invert = TRUE)
    useful_df=ldply(opt[useful_opts])
    useful_df$rename=str_replace_all(useful_df$.id, pattern = '_', replacement = '.')
    useful_df$V1=str_replace_all(useful_df$V1, pattern = '_', replacement = '.')
    useful_df[useful_df$V1 == 'FALSE','rename'] = str_replace_all(useful_df[useful_df$V1 == 'FALSE','rename'], 'add', 'rm')
    other_opts = grep('TRUE|FALSE', useful_df$V1, invert = TRUE)
    useful_df[other_opts, 'rename'] = paste0(useful_df[other_opts, 'rename'], '.', useful_df[other_opts, 'V1'])
    #str(useful_df)
    dir_name = paste(useful_df$rename, collapse = '_')
    return (dir_name)
}


#input_data = mean_prediction
f_plot_model_data <- function(input_data, plot_title, save_path){
    cat('In plot model', plot_title, '\n')
    input_data$residual = mean_prediction$obs - mean_prediction$pred
    #str(input_data)
    library(gridExtra)
    hist_plot <-ggplot(input_data, aes(obs)) + geom_histogram(binwidth = 0.2, colour = "black", fill = "white") + facet_wrap(~population)
    density_plot <- ggplot(input_data, aes(obs, color = population)) + geom_histogram(binwidth = 0.2, colour = "black", fill = "white") + ggtitle(plot_title)
    x_range = range(input_data$obs)
    pred_vs_obs <- ggplot(input_data, aes(x = pred, y = obs, color = population)) + geom_point() + theme(legend.position="top") + xlim(range(input_data$obs))
    residual_plot<-ggplot(input_data, aes(x = pred, y = residual, color = population)) + geom_point() + theme(legend.position="top")
    final_plot <-arrangeGrob(density_plot, hist_plot ,pred_vs_obs, residual_plot, ncol=2, nrow =2)
    #ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
    #ggsave(filename = save_path, final_plot)

    tiff(save_path)
    plot(final_plot)
    dev.off()
    #class(final_plot)
    return (final_plot)
    #plot(final_plot)

}


f_get_TF_expression<- function(output_dir, type = 'TF'){
    if(type=='TF'){
        tf_gene_expression = read.table(f_p('%s/rnaseq/tf_transcript_data.bed', output_dir), header = T, stringsAsFactors = FALSE, na.strings = 'NA')
    }else if(type == 'fakeTF'){
        tf_gene_expression = read.table(f_p('%s/rnaseq/faked_tf_transcript_data.bed', output_dir), header = T, stringsAsFactors = FALSE, na.strings = 'NA')
    }else if(type == 'miRNA'){
        tf_gene_expression = read.table(f_p('%s/rnaseq/random_miRNA_transcript_data.bed', output_dir), header = T, stringsAsFactors = FALSE, na.strings = 'NA')
    }

    return (tf_gene_expression)
}



f_add_tf_interactions <- function(data){
                                        #Add the TF interactions
    tf_regions = data
    tf_regions = tf_regions[grep('DNase|H[0-9]K[0-9]|RNASEQ', tf_regions$feature, invert= TRUE),]
    tf_regions=tf_regions[!duplicated(tf_regions),]

    genome_ragnes = makeGRangesFromDataFrame(tf_regions[,c('chr', 'feature_start','feature_end')])

    matches = as.data.frame( findOverlaps(genome_ragnes, genome_ragnes, minoverlap = 200) )

    matches = matches[matches$queryHits != matches$subjectHits, ]
    
    if(nrow(matches) > 0){
                                        #f_one_pair_tf_interaction(match_line, sample_cols, tf_regions)    
        
        overlap_df=data.frame(matches)
        str(overlap_df)
        head(overlap_df)
        overlap_df$name = paste0(tf_regions[overlap_df$queryHits,'feature_tf'], '-', tf_regions[overlap_df$subjectHits,'feature_tf'] )
        overlap_pairs = overlap_df[overlap_df$name %in% valid_interaction$V1,]
        flog.info('%s out of %s is valiad TF interactions',nrow(overlap_pairs), nrow(matches))
                                        #str(overlap_df)
        if(nrow(overlap_pairs)>0){
            tf_interaction_impact = ldply(  apply(overlap_pairs[,1:2], MARGIN = 1, f_one_pair_tf_interaction, sample_cols, tf_regions) )
            dim(tf_interaction_impact)
            head(tf_interaction_impact)
            
            tf_interaction_impact$.id = NULL
            row.names(tf_interaction_impact) = paste0('TF.overlap.',make.names(tf_interaction_impact$feature, unique = TRUE))
            tf_valid_interaction_impact = tf_interaction_impact[tf_interaction_impact$feature %in% valid_interaction$V1,]
            cat('Interaction terms', dim(tf_valid_interaction_impact), '\n')
            
            transcript_data_merge = rbind(data, tf_valid_interaction_impact)
            
        }else{
            cat('Empty overlaps','\n')
        }
        
    }else{
        transcript_data_merge = data
    }

    
   ############
   #Promoter - enhancer TF interactions
   ############

    promoter_index=as.numeric(which(tf_regions$type == 'promoter'))
    enhancer_index=as.numeric(which(tf_regions$type == 'enhancer'))
    pair_df=data.frame(promoter = rep(promoter_index, each = length(enhancer_index)), enhancer = rep(enhancer_index, times = length(promoter_index)))
                                        #str(pair_df)
    if(nrow(pair_df) > 0){
        pair_df$name = paste0(tf_regions[pair_df$promoter,'feature_tf'], '-', tf_regions[pair_df$enhancer,'feature_tf'] )
        promoter_pairs=pair_df[pair_df$name %in% valid_interaction$V1,]
        str(promoter_pairs)

        if(nrow(promoter_pairs) > 0){
            
            promoter_interaction_impact = ldply(  apply(promoter_pairs[,1:2], MARGIN = 1, f_one_pair_tf_interaction, sample_cols, tf_regions) )
            dim(promoter_interaction_impact)
            promoter_interaction_impact$.id=NULL
            rownames(promoter_interaction_impact) = make.names(paste0('P_E.', promoter_interaction_impact$feature), unique = TRUE)
            
            dim(transcript_data_merge)
            transcript_data_merge = rbind(transcript_data_merge, promoter_interaction_impact)
            cat('Interaction terms', dim(promoter_interaction_impact), '\n')
            rm(promoter_interaction_impact)
        }else{
            cat('Empty promoter-enhancers!', '\n')
        }
    }

    return (transcript_data_merge)
    #data <<- transcript_data_merge
}


f_add_related_miRNAs <- function(transcript_id, input_data){
    miRNA_target_table = read.table('./data/raw_data/miRNA/miRNA_ensemble.txt', header = TRUE)
    miRNA_expression = read.table('./data/raw_data/miRNA/GD452.MirnaQuantCount.1.2N.50FN.samplename.resk10.txt', header = TRUE)
    related_miRNAs = subset(miRNA_target_table, ensembl_gene_id == str_split(transcript_id, '[.]')[[1]][1] )$miRNA
    cat('Transcript', transcript_id, '\n')
    miRNA_expression$TargetID = as.character(miRNA_expression$TargetID)
    correlated_miRNA_expression = subset(miRNA_expression, TargetID %in% related_miRNAs )[, c('TargetID', sample_cols)]
    row.names(correlated_miRNA_expression) = correlated_miRNA_expression$TargetID
    cat('Add miRNA to the training:', as.character(correlated_miRNA_expression$TargetID), '\n')
    correlated_miRNA_expression$TargetID = NULL
    output_data = cbind(input_data, t(correlated_miRNA_expression[rownames(final_train_data)]))
    return (output_data)
}



f_add_predicted_TF_expression <- function(add_predict_TF, batch_name, input_data){

    output_data = input_data
    if (add_predict_TF == TRUE){

        cat('======== add predict TF ========\n')
        
        tf_prediction = f_p('./data/%s/rnaseq/chrTF/tf_prediction', batch_name)

        tf_gene_id_unique=tf_gene_id[!duplicated(tf_gene_id$ensembl_gene_id),]
        rownames(tf_gene_id_unique) = tf_gene_id_unique$ensembl_gene_id

        
        tf_prediction_data = read.table(tf_prediction, header = TRUE, na.strings = 'NA')
        
        rownames(tf_prediction_data) = paste0( 'tf_predict', tf_gene_id_unique[tf_prediction_data$gene,'tf'])
        
        output_data = cbind(input_data, data.frame(t(tf_prediction_data[, rownames(input_data) ])))
        
    }

    return (output_data)

}





f_get_test_data <- function(empty_data = TRUE){
    fun_dir = './data/test_data/s_gene_regression_fun'
    
    #gene = 'ENSG00000241973.6'
    #batch_name = '462samples_quantile_rmNA'

    #Before add the TF concentration
    #write.table(transcript_data, f_p('%s/transcript_data', fun_dir))
    if (empty_data == FALSE){
        test_data = read.table(f_p('%s/transcript_data', fun_dir))
        head(test_data)
    }else{
        test_data = NULL
    }

    sample_info = read.table(f_p('./data/462samples/chr_vcf_files/integrated_call_samples_v3.20130502.ALL.panel'), header = TRUE, row.names = 1)
    snyder_samples = read.table(f_p('./data/54samples_evalue/output/sample.list'),
                           header = FALSE)
    
    return (list(data=test_data, sample_info = sample_info, snyder_samples = snyder_samples$V2))
}




f_correct_gm12878_bias <- function(input_data, sample_id, exclude_rows = 'gene.RNASEQ'){

    sample_cols = sort(grep('(NA|HG)[0-9]+', colnames(input_data),value = T))

    #head10(input_data)
    correct_data = input_data

    selected_rows = grep(exclude_rows, rownames(input_data), invert = TRUE, value = TRUE)

    
    if (sample_id %in% colnames(input_data)){
        correct_data[selected_rows, sample_cols] = input_data[selected_rows, sample_cols] - input_data[selected_rows, sample_id]
    }else{
        flog.error('%s is not in the colonames of input data', sample_id)
        return (input_data)
    }
    
    return (correct_data)
}



t_correct_gm12878_bias <-function(){
    input_data = f_get_test_data()
    sample_id = 'NA12872'
    
    corrected_data = f_correct_gm12878_bias(input_data$data, sample_id)

    f_assert(all(corrected_data[-1, sample_id] == 0), 'correction is wrong')
}



f_add_population_and_gender <- function(obj, add_YRI, select_pop){
    #input_data$population = NULL
    
    input_data = obj$data
    snyder_samples = obj$snyder_samples
    input_data[,'population'] = as.character(sample_info[rownames(input_data), c('pop')])
    input_data[,'gender'] = sample_info[rownames(input_data), c('gender')] == 'male'
    input_data$population = as.factor(input_data$population)
    input_data$gender = as.factor(input_data$gender)
 
    #Remove the population
    if (add_YRI == FALSE){
        flog.info('Remove the YRI')
        output_data = subset(input_data, population != 'YRI')
    }

    shared_cols = intersect(snyder_samples, rownames(input_data))
    
    if(select_pop == 'YRI'){
        output_data = subset(input_data, population == select_pop)
        flog.info('Subset to %s size: %s', select_pop, dim(input_data))
    }else if(select_pop == '54snyder'){
        output_data = input_data[shared_cols, ]
    }
    else if(select_pop == '29snyder'){
        
        output_data = input_data[sample(shared_cols, size = 29), ]
        
    }else if(select_pop == '29YRI'){
        YRI_samples = rownames( subset(obj$sample_info, pop == 'YRI') )
        length(YRI_samples)

        other_YRI_samples = setdiff( intersect(rownames(input_data), YRI_samples), shared_cols)
        output_data = input_data[other_YRI_samples,]
        
    }else if(select_pop == 'all'){
        output_data = input_data
    }else if(select_pop == 'None'){
        output_data = input_data
        output_data$population = NULL
        output_data$gender = NULL
    }
    else{
        flog.warn('Unknown population %s', select_pop)
    }

    return (output_data)
}




 t_add_population_and_gender <-function(){
    obj = f_get_test_data()
    obj$data = as.data.frame(t(obj$data))
    snyder_data = f_add_population_and_gender(obj, add_YRI = TRUE, select_pop = '54snyder')
    f_assert(nrow(snyder_data)== 53, 'Test sub-population')
    dim(snyder_data)
    

    YRI29_data = f_add_population_and_gender(obj, add_YRI = TRUE, select_pop = '29YRI')
    f_assert(nrow(YRI29_data)== 29, 'Test sub-population')


    snyder29_data = f_add_population_and_gender(obj, add_YRI = TRUE, select_pop = '29snyder')
    f_assert(nrow(snyder29_data)== 29, 'Test sub-population')
    
    f_assert( length(intersect(rownames(YRI29_data), rownames(snyder_data) )) == 0, 'YRI selection' )
    
}



library(RUnit)
options(run.main = F)
if (getOption('run.main', default=TRUE)) {
  runTestFile('./s_gene_regression_fun.R', testFuncRegexp = '^t_.*')
}









