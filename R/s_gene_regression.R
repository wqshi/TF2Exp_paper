
#install.packages('bnlearn')
library(bnlearn)
source('~/R/s_function.R', chdir = TRUE)
library(doMC)

library(stringr)
#install.packages('doMC')


#my_train = final_train_data
#target_col = 'gene.RNASEQ'


f_caret_regression_return_fit <- function(my_train, target_col, learner_name, quite=FALSE)
{
  
  #cl <- makeCluster(3)
  #registerDoSNOW(cl)
    registerDoMC(cores = 4)
    library(caret)
    
    dim(my_train)
    head(my_train)
    zero_cols = nearZeroVar(my_train, freqCut = 90/10)
    if (length(zero_cols) == 0)
        non_zero_data = my_train
    else
        non_zero_data = my_train[, -zero_cols]
    #non_zero_data = my_train
    cat('After filtering, training size', dim(non_zero_data), '(features, samples)')
    colnames(non_zero_data)
    
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = 5,
                             repeats = 1)
  
  if (learner_name == 'enet')
  {
    enetGrid <- expand.grid(.lambda = c(0, 0.01, .1), .fraction = seq(.05, 1, length = 5 ))

    
    Fit <- train(form = as.formula(f_p('%s ~ .', target_col)), data = non_zero_data,
                 method = "enet",
                 #tuneGrid = enetGrid,
                 trControl = fitControl,
                 preProc = c("center", "scale"))
    if (quite == FALSE){
        print(Fit$resample)
        summary(Fit)
    }


    ###Get the coefficient from the model########
    
    coeff=predict(Fit$finalModel, type='coefficient', s=Fit$bestTune$fraction, mode='f')
    print(coeff$coefficients[coeff$coefficients != 0])
    
    #train_pred = predict(Fit, newdata = my_train, type='raw')
    
    #processed=preProcess(my_train, method = c('scale','center'))
    #processed_data=predict(processed, my_train)
    #result= data.frame(truth  = my_train[[target_col]])
    #result[,predictors(Fit)] =   processed_data[,predictors(Fit)]
    
  }else
  {
    Fit <- caret::train( as.formula(f_p('%s ~ .', target_col)), data = non_zero_data, metric = "Rsquared", 
                         method = learner_name,
                         trControl = fitControl,
                         verbose = TRUE)
    
    #rfFit = randomForest(as.formula(f_p('%s ~ .', target_col)), data = my_train, ntree = 1000, do.trace=100)
    
    #train_pred = rfFit$pred
    #result= data.frame(truth  = my_train[[target_col]])
  }
  
  
  
  return (Fit)
  
}



#install.packages('optparse')
library("optparse")

option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="dataset file name", metavar="character"),
    make_option(c("--add_histone"), type="character", default='TRUE', help="output file name [default= %default]", metavar="character"),
    make_option(c("--add_miRNA"), type="character", default='TRUE', help="Add miRNA or not", metavar="character"),
    make_option(c("--test"), type="character", default=NULL, help="output file name [default= %default]", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$batch_name)){
    batch_name = '54samples_peer'
    add_histone = TRUE
    add_miRNA = TRUE
    test_flag = TRUE
}else{
    batch_name = opt$batch_name
    add_histone = opt$add_histone == 'TRUE'
    add_miRNA = opt$add_miRNA == 'TRUE'
    test_flag = opt$test == 'TRUE'
    cat('Test flag :', test_flag, '\n')
}

chr_str = 'chr22'
model_str = 'all'


output_dir = f_p('./data/%s/', batch_name)
expression_file = paste0(output_dir, 'rnaseq/gene_regulatory_region_feature_profile.', chr_str )

expression_data = read.table(expression_file, header = TRUE, na.strings = '.')
head(expression_data)
tail(expression_data[,1:20])



cat('Number regulatory regions and empty fragment ids', '\n')
table(expression_data$type, expression_data$hic_fragment_id == '.')

table(expression_data$type, expression_data$feature)

table(unlist(table(expression_data$gene)))

colnames(expression_data)


#Get the expressed transcripts
sample_cols = grep('(NA|HG)[0-9]+', colnames(expression_data),value = T)


expression_data = expression_data[,c('chr', 'start', 'end', 'gene','feature','feature_start', 'feature_end', 'hic_fragment_id' ,'type', sample_cols )]
#write.table(sample_cols, './data/output/sample.list', sep = '\t', quote = F, col.names = F)
rna_data = subset(expression_data, feature == 'RNASEQ')
rownames(rna_data) = rna_data$gene
expressed_transcripts = rownames(rna_data)[rowSums(rna_data[, sample_cols]) > 0]
cat('The number of expressed transcripts:', length(expressed_transcripts), '\n')

head(rna_data)
#str(rna_data)
#Get the genes with at least one elements
gene_with_regulatory_elements = as.data.frame.matrix(table(expression_data$gene, expression_data$feature))
head(gene_with_regulatory_elements)
str(gene_with_regulatory_elements, max.level = 1)
regulated_genes_names = rownames(gene_with_regulatory_elements)[ rowSums(gene_with_regulatory_elements) > 10]
cat('The number of regulated transcripts:', length(regulated_genes_names), '\n')

#Get the inverstigeed genes
genes_names = intersect(expressed_transcripts, regulated_genes_names)
cat('The number of investigated transcripts:', length(genes_names), '\n')


length(genes_names)

print(genes_names[1:20])

prediction_performance = data.frame(gene = 'mean', performance = '0', stringsAsFactors = F)

sample_info = read.table(f_p('%s/chr_vcf_files/integrated_call_samples_v3.20130502.ALL.panel', output_dir ), header = TRUE, row.names = 1)


sample_cols = grep('(NA|HG)[0-9]+', colnames(expression_data), value = TRUE)

#Read the miRNA interaction and expression data.

miRNA_target_table = read.table('/homed/home/shi/expression_var/data/raw_data/miRNA/miRNA_ensemble.txt', header = TRUE)
miRNA_expression = read.table('/homed/home/shi/expression_var/data/raw_data/miRNA/GD452.MirnaQuantCount.1.2N.50FN.samplename.resk10.txt', header = TRUE)

i = 1

for (i in 1:length(genes_names)){
    if(i == 50 & test_flag == TRUE){
        print( 'First 500 genes')
        break
    }
    
    transcript_id = genes_names[i]
    cat('\n','============',i, transcript_id, '==============','\n')
    sum(gene_with_regulatory_elements[transcript_id,])
    gene_with_regulatory_elements[transcript_id,]
    transcript_data = expression_data[expression_data$gene == transcript_id,]
    dim(transcript_data)

    
    head(transcript_data)
    table(transcript_data$feature)
    train_data = transcript_data[,sample_cols]
    sum(is.na(train_data))
    dim(train_data)

    head(train_data)
    
    train_data[is.na(train_data)] = 0
    rownames(train_data) = make.names(paste0(transcript_data$type, '-',transcript_data$feature), unique = T)
    train_data_rmdup = train_data[!duplicated(train_data),]

    train_data2 = as.data.frame(t(train_data_rmdup))
    head(train_data2)
    if (add_histone == FALSE){
        none_histone_cols = grep('H3K',colnames(train_data2), value = TRUE, invert = TRUE)
        final_train_data = train_data2[, none_histone_cols]
    }else{
        final_train_data  = train_data2
    }


    head(miRNA_target_table)
    dim(train_data2)
    
    related_miRNAs = subset(miRNA_target_table, ensembl_gene_id == str_split(transcript_id, '[.]')[[1]][1] )$miRNA

    if(length(related_miRNAs) > 0){
        cat('Transcript', transcript_id, '\n')
        miRNA_expression$TargetID = as.character(miRNA_expression$TargetID)
        correlated_miRNA_expression = subset(miRNA_expression, TargetID %in% related_miRNAs )[, c('TargetID', sample_cols)]
        row.names(correlated_miRNA_expression) = correlated_miRNA_expression$TargetID
        cat('Add miRNA to the training:', as.character(correlated_miRNA_expression$TargetID), '\n')
        correlated_miRNA_expression$TargetID = NULL
        final_train_data = cbind(final_train_data, t(correlated_miRNA_expression[rownames(final_train_data)]))
    }

    #ead(final_train_data)
    colnames(final_train_data)
    
    head(final_train_data)
    rownames(final_train_data)
    final_train_data$population = NULL
    final_train_data[,'gender'] = sample_info[rownames(final_train_data), c('gender')] == 'male'
    #final_train_data$population = as.factor(final_train_data$population)
    #final_train_data$gender = as.factor(final_train_data$gender)
    table(final_train_data$population)
    
    result <- try(
                        fit  <- f_caret_regression_return_fit(my_train = final_train_data, target_col = 'gene.RNASEQ', learner_name = 'enet'),
                        silent=TRUE
                       )
    if (class(result) == "try-error"){
        cat('Error in ', transcript_id, '\n')
        next
    }

    
    #fit  <- f_caret_regression_return_fit(my_train = final_train_data, target_col = 'gene.RNASEQ', learner_name = 'glmnet')

    head(final_train_data)
    
    #print(fit$results)

    max_performance = max(fit$results$Rsquared)
    print(max_performance)
    prediction_performance = rbind(prediction_performance, c(as.character(transcript_id), as.character(max_performance)))
    #print(prediction_performance)
}


output_file = f_p('%s/output/lasso.prediction.%s.%s.%s.%s.txt', output_dir, chr_str, model_str, ifelse(add_histone,  'addHis', 'rmvHis'), ifelse(add_miRNA, 'addMiR', 'rmvMiR') )
prediction_performance= prediction_performance[-1,]
prediction_performance$performance = as.numeric(prediction_performance$performance)
prediction_performance['mean',] =c('mean', mean(prediction_performance[complete.cases(prediction_performance),]$performance))

write.table(prediction_performance, file = output_file , sep = '\t', quote=FALSE, row.names = TRUE)

performance_data = read.table(file = output_file, sep = '\t')
cat("Mean R-square:",  mean(performance_data[complete.cases(performance_data),]$performance), '\n')












