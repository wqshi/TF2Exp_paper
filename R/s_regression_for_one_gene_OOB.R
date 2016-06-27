##This version is try to used OOB pro

setwd('~/expression_var/R/')
source('s_gene_regression_fun.R')
library(stringr)
source('s_project_funcs.R')
source('s_gene_data_class.R')
library(futile.logger)
library("optparse")
library(plyr)
library(dplyr)
flog.info('In OOB')

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
    make_option(c("--add_gm12878"), type="character", default='TRUE', help="Whether to remove the GM12878 from the TF impact matrix", metavar="character"),
    make_option(c("--new_batch"), type="character", default='', help="Change the expression data to another batch", metavar="character"),
    make_option(c("--batch_mode"), type="character", default='TF', help="Change the expression data to another batch", metavar="character")
);

flog.info('Before the opt parse')
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

str(opt)

if (is.null(opt$batch_name)){
    batch_name = '54samples_evalue'
    #batch_name = '462samples_genebody'
    #batch_name = '462samples_quantile_rmNA'
    #batch_name = 'test'
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
    #gene = 'ENSG00000093072.11' # Gene in the test dir.


    #gene = 'ENSG00000100321.10' #Gene with big difference when remove self interaction
    #gene = 'ENSG00000184702.13' #with smaler size.
    permutation_flag = FALSE
    chr_str = 'chr22'
    add_TF_exp_only=FALSE
    add_predict_TF=FALSE
    add_YRI=TRUE
    select_pop= 'all' #'None' #'all'
    TF_exp_type = 'fakeTF'
    add_gm12878=FALSE
    new_batch = '54samples_peer'
    #new_batch = '462samples_genebody'
    batch_mode = 'TF'
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
    new_batch = opt$new_batch
    batch_mode = opt$batch_mode 
    tuneGrid = tuneGridList[[train_model]]
}

flog.info('Begin')

model_str = 'all'
cat('Tune Grid is NULL:', is.null(tuneGrid), '\n')
print(tuneGrid)


output_dir = f_p('./data/%s/', batch_name)

test_gene <- GENE(data = data.frame(), gene_name = gene, chr_str = chr_str, batch_name = batch_name)
test_gene$read_data()
sample_cols = test_gene$get_samples()
test_gene$change_expression(new_batch, batch_mode)
test_gene$subset_snps_in_tf_regions(batch_mode)

expression_data = test_gene$data

dim(expression_data)
colnames(expression_data)

table(expression_data$type)

gene_with_regulatory_elements = as.data.frame.matrix(table(expression_data$gene, expression_data$feature))
#head(gene_with_regulatory_elements)
regulated_genes_names = rownames(gene_with_regulatory_elements)[ rowSums(gene_with_regulatory_elements) > 10]
cat('The number of regulated transcripts:', length(regulated_genes_names), '\n')

#Get the inverstigeed genes
genes_names = regulated_genes_names
cat('The number of investigated transcripts:', length(genes_names), '\n')

prediction_performance = data.frame(gene = 'mean', train_performance=0, performance = '0', SD = '0', num_feature = '0' ,stringsAsFactors = F)
collected_features = data.frame()
sample_info = read.table(f_p('%s/chr_vcf_files/integrated_call_samples_v3.20130502.ALL.panel', output_dir ), header = TRUE, row.names = 1)

#Read the miRNA interaction and expression data.
#miRNA_expression = read.table('./data/raw_data/miRNA/GD452.MirnaQuantCount.1.2N.50FN.samplename.resk10.txt', header = TRUE)


non_sample_cols = setdiff(colnames(expression_data), sample_cols)
#sample_cols = intersect(sample_cols, colnames(miRNA_expression))


#Read the TF expression data.
tf_gene_id = read.table('./data/raw_data/rnaseq/tf_ensemble_id.txt', header = T, stringsAsFactors = FALSE)
rownames(tf_gene_id) = tf_gene_id$tf
expression_data$feature_tf =''
expression_data$feature_tf = (tf_gene_id[as.character(expression_data$feature),'external_name'])

tf_gene_expression = f_get_TF_expression(output_dir, type = TF_exp_type)
rownames(tf_gene_expression) = tf_gene_id$tf
flog.info('After load the data')

opt_name = f_convet_opts_to_output_dir(opt)
results_dir = f_p('%s/rnaseq/%s/%s/', output_dir, chr_str, opt_name )
dir.create(results_dir, showWarnings = FALSE)
flog.info('Results dir: %s', opt_name)
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
    train_data_raw = f_add_tf_interactions_old2(transcript_data, debug = F)
    train_data = train_data_raw[,sample_cols]
    train_data[is.na(train_data)] = 0
    train_data_rmdup = train_data[!duplicated(train_data),]
    rownames(train_data)
    
    
    train_data2 = as.data.frame(t(train_data_rmdup))
    colnames(train_data2)
    if (add_histone == FALSE){
        none_histone_cols = grep('H3K',colnames(train_data2), value = TRUE, invert = TRUE)
        final_train_data = train_data2[, none_histone_cols]
    }else{
        final_train_data  = train_data2
    }    
    dim(final_train_data)
    colnames(final_train_data)
    if(add_miRNA == TRUE){
        final_train_data = f_add_related_miRNAs(transcript_id, final_train_data)
    }
    ##Make a TF concentration only model.
    tf_expression = data.frame(t(tmp[, rownames(final_train_data) ]))
    tf_expression = tf_expression[,!duplicated(t(tf_expression))]
    #str(tf_expression)
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

    
    dim(final_train_data)
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
    
    final_train_data <- f_filter_training_features(final_train_data, batch_mode)    
    flog.info('After add other feature data')

    #penalty_factors = f_add_penalty_factor(final_train_data, train_data_raw[,c('chr','cor')])
    #dim(final_train_data)
    #length(penalty_factors)
    
    output_figure_path = f_p('%s/%s.learn.curve.tiff',  results_dir, gene)
    result <- try(
        fit  <- f_caret_regression_return_fit(my_train = final_train_data, target_col = 'gene.RNASEQ', learner_name = train_model, tuneGrid, output_figure_path)
       ,silent=TRUE)
        
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

    train_perf = f_get_train_performance(fit, final_train_data, 'gene.RNASEQ')
    f_plot_model_data(mean_prediction, plot_title = f_p('%s \n R-squared: %.2f', transcript_id, fit$max_performance), save_path = f_p('%s/%s.enet.tif', results_dir, transcript_id) )

    #plot(tmp_plot)
    #ggsave('a.tiff', tmp_plot)
    #ggsave('a2.tiff', tmp_plot, dpi = 72)
    
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

    collected_features = rbind(collected_features, features_df)
    
 
    max_performance = fit$max_performance
    RsquaredSD =fit$results[rownames(fit$bestTune), 'RsquaredSD' ]
    print(fit$results)
    cat("\n", 'train_perf', train_perf ,'performance', max_performance, 'SD', RsquaredSD, 'Number of features:', nrow(features_df) )
    prediction_performance = rbind(prediction_performance, c(as.character(transcript_id), as.character(train_perf), as.character(max_performance), as.character(RsquaredSD), as.character(fit$num_features)))
    #print(prediction_performance)
}



output_file = f_p('%s/%s.enet',  results_dir, gene)

flog.info('Output file: %s', output_file)

prediction_performance= prediction_performance[-1,]
write.table(prediction_performance, file = output_file , sep = '\t', quote=FALSE, row.names = FALSE)
write.table(collected_features, file = f_p('%s.features', output_file), sep = '\t', quote = FALSE, row.names =FALSE)
