setwd('~/expression_var/R/')
source('s_gene_regression_fun.R')
source('s_project_funcs.R')
library(stringr)

library(futile.logger)
library("optparse")
library(plyr)
library(dplyr)


option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="dataset file name", metavar="character"),
    make_option(c("--add_histone"), type="character", default='TRUE', help="output file name [default= %default]", metavar="character"),
    make_option(c("--add_miRNA"), type="character", default='FALSE', help="Add miRNA or not", metavar="character"),
    make_option(c("--add_TF_exp"), type="character", default='FALSE', help="Add TF expression levels into the model or not", metavar="character"),
    make_option(c("--test"), type="character", default=NULL, help="output file name [default= %default]", metavar="character"),
    make_option(c("--gene"), type="character", default='', help="The name of gene", metavar="character"),
    make_option(c("--model"), type="character", default='enet', help="The machine learning method used, eg. enet, and rfe", metavar="character"),
    make_option(c("--add_permutation"), type="character", default='FALSE', help="Permutate the train features", metavar="character"),
    make_option(c("--chr_str"), type="character", default='chr22', help="Chromosome name", metavar="character"),
    make_option(c("--add_TF_exp_only"), type="character", default='FALSE', help="Add TF expression as the input features", metavar="character"),
    make_option(c("--add_predict_TF"), type="character", default='FALSE', help="Add TF predicted expression as the input features", metavar="character"),
    make_option(c("--add_YRI"), type="character", default='FALSE', help="Whether to remove YRI population ", metavar="character"),
    make_option(c("--population"), type="character", default='all', help="Whether to remove YRI population ", metavar="character"),
    make_option(c("--TF_exp_type"), type="character", default='TF', help="Read TF expression, faked TF, or random miRNA", metavar="character"),
    make_option(c("--add_gm12878"), type="character", default='TRUE', help="Whether to remove the GM12878 from the TF impact matrix", metavar="character")
);

flog.info('Before the opt parse')
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

str(opt)


if (is.null(opt$batch_name)){
    batch_name = '54samples_evalue'
    #batch_name = '462samples_sailfish_quantile'
    #batch_name = '462samples_quantile_rmNA'
    add_histone = TRUE
    add_miRNA = FALSE
    add_TF_exp = FALSE
    test_flag = TRUE
    tuneGrid = NULL
    #train_model = 'rf'
    train_model = 'glmnet'
    #gene = 'ENSG00000235478.1'
    gene = 'ENSG00000241973.6'# The gene with old recode and good accruacy
    #gene='ENSG00000196576.10' : largest memory
    permutation_flag = FALSE
    chr_str = 'chr22'
    add_TF_exp_only=TRUE
    add_predict_TF=FALSE
    add_YRI=TRUE
    select_pop='all'
    TF_exp_type = 'fakeTF'
    add_gm12878=FALSE
}else{
    batch_name = opt$batch_name
    add_histone = opt$add_histone == 'TRUE'
    add_miRNA = opt$add_miRNA == 'TRUE'
    add_TF_exp = opt$add_TF_exp == 'TRUE'
    permutation_flag = opt$add_permutation == 'TRUE'
    test_flag = opt$test == 'TRUE'
    train_model = opt$model
    gene=opt$gene
    output_mode = opt$output_mode
    cat('Test flag :', test_flag, 'Model ', train_model, '\n')
    chr_str = opt$chr_str
    add_TF_exp_only=opt$add_TF_exp_only == 'TRUE'
    add_predict_TF= opt$add_predict_TF == 'TRUE'
    add_YRI = opt$add_YRI == 'TRUE'
    select_pop=opt$population
    TF_exp_type=opt$TF_exp_type
    add_gm12878=opt$add_gm12878 == 'TRUE'
    tuneGrid = tuneGridList[[train_model]]
}

flog.info('Begin')

model_str = 'all'
cat('Tune Grid is NULL:', is.null(tuneGrid), '\n')
print(tuneGrid)


output_dir = f_p('./data/%s/', batch_name)
expression_file = f_p('%s/rnaseq/%s/%s.txt', output_dir, chr_str, gene)

expression_data = read.table(expression_file, header = TRUE, na.strings = 'NA')
head(expression_data[,1:10])
#tail(expression_data[,1:20])


cat('Number regulatory regions and empty fragment ids', '\n')
#table(expression_data$type, expression_data$hic_fragment_id == '.')

#table(expression_data$type, expression_data$feature)
#colnames(expression_data)


#Get the expressed transcripts
sample_cols = sort(grep('(NA|HG)[0-9]+', colnames(expression_data),value = T))
cat('The number of samples:', length(sample_cols), '\n')
expression_data = expression_data[,c('chr', 'start', 'end', 'gene','feature','feature_start', 'feature_end', 'hic_fragment_id' ,'type', sample_cols )]
#write.table(sample_cols, './data/output/sample.list', sep = '\t', quote = F, col.names = F)
head10(expression_data)
gene_with_regulatory_elements = as.data.frame.matrix(table(expression_data$gene, expression_data$feature))
head(gene_with_regulatory_elements)
regulated_genes_names = rownames(gene_with_regulatory_elements)[ rowSums(gene_with_regulatory_elements) > 10]
cat('The number of regulated transcripts:', length(regulated_genes_names), '\n')

#Get the inverstigeed genes
genes_names = regulated_genes_names
cat('The number of investigated transcripts:', length(genes_names), '\n')

prediction_performance = data.frame(gene = 'mean', performance = '0', SD = '0', num_feature = '0' ,stringsAsFactors = F)
collected_features = data.frame()
sample_info = read.table(f_p('%s/chr_vcf_files/integrated_call_samples_v3.20130502.ALL.panel', output_dir ), header = TRUE, row.names = 1)
if ( length(grep('random', batch_name)) != 0  ){
    rownames(sample_info) = sample(rownames(sample_info), size = nrow(sample_info))
}

#Read the miRNA interaction and expression data.
miRNA_target_table = read.table('./data/raw_data/miRNA/miRNA_ensemble.txt', header = TRUE)
miRNA_expression = read.table('./data/raw_data/miRNA/GD452.MirnaQuantCount.1.2N.50FN.samplename.resk10.txt', header = TRUE)


non_sample_cols = setdiff(colnames(expression_data), sample_cols)
sample_cols = intersect(sample_cols, colnames(miRNA_expression))


#Read the TF expression data.
tf_gene_id = read.table('./data/raw_data/rnaseq/tf_ensemble_id.txt', header = T, stringsAsFactors = FALSE)
rownames(tf_gene_id) = tf_gene_id$tf
expression_data$feature_tf =''
expression_data$feature_tf = (tf_gene_id[expression_data$feature,'external_name'])

tf_gene_expression = f_get_TF_expression(output_dir, type = TF_exp_type)
rownames(tf_gene_expression) = tf_gene_id$tf

valid_interaction = read.table('./data/raw_data/biogrid/tf_interactions.txt')

flog.info('After load the data')

opt_name = f_convet_opts_to_output_dir(opt)
results_dir = f_p('%s/rnaseq/%s/%s/', output_dir, chr_str, opt_name )
dir.create(results_dir, showWarnings = FALSE)

for (i in 1:length(genes_names)){
    
    transcript_id = genes_names[i]

    if (str_replace(transcript_id, pattern = '[.][0-9]+', replacement = '') %in% tf_gene_id$ensembl_gene_id){
        cat('Skip TF gene', transcript_id, '\n')
    }

    
    cat('\n','============',i, transcript_id, '==============','\n')
    sum(gene_with_regulatory_elements[transcript_id,])
    gene_with_regulatory_elements[transcript_id,]
    transcript_data = expression_data[expression_data$gene == transcript_id,]
    transcript_data[1,]
    dim(transcript_data)
    head(transcript_data[,1:15])
    transcript_data = transcript_data[!duplicated(transcript_data),]
    dim(transcript_data)
    table(transcript_data$feature)
    rownames(transcript_data) = make.names(paste0(transcript_data$type, '-',transcript_data$feature), unique = T)
    transcript_data[is.na(transcript_data)]=0 #Some TF regions don't have variations 
    sum(is.na(transcript_data))
    head10(transcript_data)


    if(add_gm12878 == FALSE){
        transcript_data = f_correct_gm12878_bias(transcript_data, 'NA12872')
        flog.info('Remove NA12872 bias.')
        head10(transcript_data)
    }
    
    
    #######Add the TF concentration data#########
    tmp=tf_gene_expression[as.character(transcript_data$feature),sample_cols]
    tmp[is.na(tmp)]=1 #Set the non-TF rows to 1
    scaled_tmp=t(apply(tmp, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))) + 1)
    scaled_tmp[is.na(scaled_tmp)] = 1 #For the DNASE regions.
    dim(scaled_tmp)
    head(tf_gene_expression)
    head(scaled_tmp[,1:15])
    dim(tmp)
    
    #transcript_data = transcript_data[!duplicated(transcript_data[, non_sample_cols[1:8]]),]
    dim(transcript_data)
    head10(tmp)
    transcript_data_tf_concentration = transcript_data
    dim(transcript_data)
    dim(scaled_tmp)
    head10(transcript_data)
    
    ###Add permutation within features.#####
    if (permutation_flag == TRUE){
        sum(rowSums(is.na(transcript_data)) < 0.5 * nrow(transcript_data))
        #transcript_data[7,]        
        set.seed(11)
        features = setdiff(rownames(transcript_data), c('gene.RNASEQ'))
        permutation_index = sample(features, size = length(features))
        cat('Permuate features:', permutation_index[1:5], '\n')
        cat('Original features:', features[1:5], '\n')
        transcript_data[features, sample_cols] = transcript_data[permutation_index, sample_cols]
        
    }

    if (add_TF_exp == TRUE){
        transcript_data_tf_concentration[, sample_cols] = transcript_data[, sample_cols] * scaled_tmp[, sample_cols]
    }
    
    
    ##Add TF-TF interactions
    train_obj  = f_add_tf_interactions(transcript_data)


    
    train_data = train_obj[,sample_cols]
    train_data[is.na(train_data)] = 0
    train_data_rmdup = train_data[!duplicated(train_data),]

    
    
    train_data2 = as.data.frame(t(train_data_rmdup))
    
    if (add_histone == FALSE){
        none_histone_cols = grep('H3K',colnames(train_data2), value = TRUE, invert = TRUE)
        final_train_data = train_data2[, none_histone_cols]
    }else{
        final_train_data  = train_data2
    }    
    dim(final_train_data)

    if(add_miRNA == TRUE){
        final_train_data = f_add_related_miRNAs(transcript_id, final_train_data)
    }

    ##Make a TF concentration only model.
    tf_expression = data.frame(t(tmp[, rownames(final_train_data) ]))
    tf_expression = tf_expression[,!duplicated(t(tf_expression))]
    str(tf_expression)
    if (add_TF_exp_only == TRUE){
        cat('======== TF expression only========\n')
        tf_expression$gene.RNASEQ = final_train_data[rownames(tf_expression), 'gene.RNASEQ']
        final_train_data_bak= final_train_data
        final_train_data = tf_expression
        colnames(tf_expression)
    }else if(add_TF_exp == TRUE){
        #tf_expression = data.frame(t(tmp[, rownames(final_train_data) ]))
        cat('======== add TF expression========\n')
        final_train_data_bak= final_train_data
        #final_train_data= final_train_data_bak
        final_train_data = cbind(final_train_data, tf_expression)
    }else{
        final_train_data = final_train_data
    }

    #Add predicted TF expression
    final_train_data = f_add_predicted_TF_expression(add_predict_TF, batch_name, final_train_data)
    head10(final_train_data)
    dim(final_train_data)
    rownames(final_train_data)
    #Add population information
    obj<-f_get_test_data(empty_data = TRUE)
    obj$data = final_train_data
    final_train_data <- f_add_population_and_gender(obj, add_YRI, select_pop)
    
    additional_cols  =grep('gender|hsa-miR|population', colnames(final_train_data), value = TRUE)
    additional_data =transcript_data[rep('gene.RNASEQ', length(additional_cols)),]
    rownames(additional_data) = additional_cols
    additional_data[,sample_cols] = t(final_train_data[sample_cols, additional_cols])
    
    transcript_data = rbind(transcript_data, additional_data)
    flog.info('After add other feature data')
    #train_model = 'glmnet'
    #train_model = 'enet'
    colnames(final_train_data)
    result <- try(
                        fit  <- f_caret_regression_return_fit(my_train = final_train_data, target_col = 'gene.RNASEQ', learner_name = train_model, tuneGrid),
                        silent=TRUE
                       )

    if (class(result) == "try-error"){
        print(result)
        cat('Error in ', transcript_id, '\n')
        next
    }

    
    mean_prediction <- as.data.frame(fit$pred %>% group_by(as.factor(rowIndex)) %>% dplyr::summarise(pred = mean(pred), obs = mean(obs)))
    rownames(mean_prediction) = rownames(final_train_data)
    mean_prediction$population = final_train_data$population
    colnames(mean_prediction)
    head(mean_prediction)
    f_assert(all(final_train_data$gene.RNASEQ == mean_prediction$obs), 'Error in retrieve gene expression')
    
    grep('promoter', colnames(final_train_data), value = TRUE)
    #Write the mean prediction for each TF genes.
    if (chr_str == 'chrTF'){
        cat('load TF predictions \n')
        pred_file = f_p('%s/%s.enet.pred',  results_dir, gene)
        pred = t(mean_prediction['pred'])
        write.table(cbind(transcript_id, pred), file = pred_file, sep = '\t', quote=FALSE, row.names = FALSE)
    }

    tmp_plot<-f_plot_model_data(mean_prediction, plot_title = f_p('%s \n R-squared: %.2f', transcript_id, fit$max_performance), save_path = f_p('%s/%s.enet.tif', results_dir, transcript_id) )

    plot(tmp_plot)
    
    flog.info('After training')

    key_features = str_replace_all(names(fit$key_features),'`','')

    #transcript_data[key_features,1:10]

    feature_cols = c('chr', 'feature_start', 'feature_end')
    features_df = transcript_data[key_features, feature_cols ]
    
    rownames(features_df) = key_features
    shared_features = intersect(key_features, rownames(transcript_data))
    features_df[shared_features, ] = transcript_data[shared_features,feature_cols]

    features_df$name = paste0(transcript_id, '|' ,rownames(features_df))
    features_df$score = fit$key_features

    #print(features_df[,c('name', 'score')])

    collected_features = rbind(collected_features, features_df)
    
    #fit  <- f_caret_regression_return_fit(my_train = final_train_data, target_col = 'gene.RNASEQ', learner_name = 'glmnet')
    
    #head(final_train_data)
    #print(fit$results)
    max_performance = fit$max_performance
    RsquaredSD =fit$results[rownames(fit$bestTune), 'RsquaredSD' ]
    print(fit$results)
    cat("\n",'performance', max_performance, 'SD', RsquaredSD, 'Number of features:', nrow(features_df) )
    prediction_performance = rbind(prediction_performance, c(as.character(transcript_id), as.character(max_performance), as.character(RsquaredSD), as.character(fit$num_features)))
    #print(prediction_performance)
}



output_file = f_p('%s/%s.enet',  results_dir, gene)

flog.info('Output file: %s', output_file)

prediction_performance= prediction_performance[-1,]
write.table(prediction_performance, file = output_file , sep = '\t', quote=FALSE, row.names = FALSE)
write.table(collected_features, file = f_p('%s.features', output_file), sep = '\t', quote = FALSE, row.names =FALSE)



################
