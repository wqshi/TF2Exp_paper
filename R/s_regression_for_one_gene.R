

setwd('~/expression_var/R/')

source('~/R/s_function.R', chdir = TRUE)
source('s_gene_regression_fun.R')
library(doMC)
library(stringr)
#install.packages('doMC')
#install.packages('futile.logger')
library(futile.logger)
library(GenomicRanges)
library("optparse")
library(plyr)

option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="dataset file name", metavar="character"),
    make_option(c("--add_histone"), type="character", default='TRUE', help="output file name [default= %default]", metavar="character"),
    make_option(c("--add_miRNA"), type="character", default='FALSE', help="Add miRNA or not", metavar="character"),
    make_option(c("--add_TF_exp"), type="character", default='FALSE', help="Add TF expression levels into the model or not", metavar="character"),
    make_option(c("--test"), type="character", default=NULL, help="output file name [default= %default]", metavar="character"),
    make_option(c("--gene"), type="character", default='', help="The name of gene", metavar="character"),
    make_option(c("--model"), type="character", default='enet', help="The machine learning method used, eg. enet, and rfe", metavar="character"),
    make_option(c("--add_permutation"), type="character", default='FALSE', help="Permutate the train features", metavar="character"),
    make_option(c("--chr_str"), type="character", default='chr22', help="Chromosome name", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

str(opt)


if (is.null(opt$batch_name)){
    #batch_name = '54samples_evalue'
    batch_name = '462samples_sailfish_quantile'
    add_histone = TRUE
    add_miRNA = TRUE
    add_TF_exp = TRUE
    test_flag = TRUE
    tuneGrid = NULL
    train_model = 'glmnet'
    #gene = 'ENSG00000235478.1'
    gene = 'ENSG00000241973.6'# The gene with old recode and good accruacy
    #gene='ENSG00000196576.10' : largest memory
    permutation_flag = FALSE
    chr_str = 'chr22'
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
    tuneGrid = NULL
}

flog.info('Begin')


model_str = 'all'
cat('Tune Grid is NULL:', is.null(tuneGrid), '\n')
print(tuneGrid)


output_dir = f_p('./data/%s/', batch_name)
expression_file = f_p('%s/rnaseq/%s/%s.txt', output_dir, chr_str, gene)

expression_data = read.table(expression_file, header = TRUE, na.strings = 'NA')
#head(expression_data[,1:10])
#tail(expression_data[,1:20])


cat('Number regulatory regions and empty fragment ids', '\n')
#table(expression_data$type, expression_data$hic_fragment_id == '.')

#table(expression_data$type, expression_data$feature)
#colnames(expression_data)


#Get the expressed transcripts
sample_cols = grep('(NA|HG)[0-9]+', colnames(expression_data),value = T)

cat('The number of samples:', length(sample_cols), '\n')
expression_data = expression_data[,c('chr', 'start', 'end', 'gene','feature','feature_start', 'feature_end', 'hic_fragment_id' ,'type', sample_cols )]
#write.table(sample_cols, './data/output/sample.list', sep = '\t', quote = F, col.names = F)

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




#Read the miRNA interaction and expression data.
miRNA_target_table = read.table('./data/raw_data/miRNA/miRNA_ensemble.txt', header = TRUE)
miRNA_expression = read.table('./data/raw_data/miRNA/GD452.MirnaQuantCount.1.2N.50FN.samplename.resk10.txt', header = TRUE)

i = 1
non_sample_cols = setdiff(colnames(expression_data), sample_cols)
sample_cols = intersect(sample_cols, colnames(miRNA_expression))
length(sample_cols)

#print(non_sample_cols)

#Read the TF expression data.
tf_gene_id = read.table('./data/raw_data/rnaseq/tf_ensemble_id.txt', header = T, stringsAsFactors = FALSE)
tf_gene_id
rownames(tf_gene_id) = tf_gene_id$tf
expression_data$feature_tf =''
#str(expression_data)
head10(expression_data)


expression_data$feature_tf = (tf_gene_id[expression_data$feature,'external_name'])

table(expression_data$feature_tf)
head(expression_data$feature_tf)
table(expression_data$feature_tf)

#head(tf_gene_id)
gene_expression = read.table(f_p('%s/rnaseq/transcript_data.bed', output_dir), header = T, stringsAsFactors = FALSE, na.strings = 'NA')
rownames(gene_expression) = gene_expression$transcript_id
#str(tf_gene_id)


tf_gene_expression = gene_expression[tf_gene_id$ensembl_gene_id,]
dim(tf_gene_expression)
#head(tf_gene_expression[,1:20])

rownames(tf_gene_expression) = tf_gene_id$tf

valid_interaction = read.table('./data/raw_data/biogrid/tf_interactions.txt')

flog.info('After load the data') 

for (i in 1:length(genes_names)){
    
    transcript_id = genes_names[i]
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

    
    #######Add the TF concentration data#########
    tmp=tf_gene_expression[as.character(transcript_data$feature),sample_cols]
    tmp[is.na(tmp)]=1 #Set the non-TF rows to 1
    scaled_tmp=t(apply(tmp, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))) + 1)
    scaled_tmp[is.na(scaled_tmp)] = 1 #For the DNASE regions.

    head(tf_gene_expression)
    head(scaled_tmp[,1:15])
    dim(tmp)
    
    #transcript_data = transcript_data[!duplicated(transcript_data[, non_sample_cols[1:8]]),]
    dim(transcript_data)
    
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
        transcript_data[1,]
        transcript_data_tf_concentration[, sample_cols] = transcript_data[, sample_cols] * scaled_tmp[, sample_cols]
        #transcript_data_tf_concentration[7,]
        rowSums( transcript_data_tf_concentration == 0)
    }
    
    
    #head(transcript_data[,1:10])
    ################
    
    ##Add TF-TF interactions
    tf_regions = transcript_data_tf_concentration
    tf_regions = tf_regions[grep('DNase|H[0-9]K[0-9]|RNASEQ', tf_regions$feature, invert= TRUE),]
    tf_regions=tf_regions[!duplicated(tf_regions),]

    table(tf_regions$feature_tf)
    head(tf_regions)
    dim(tf_regions)
    
    colnames(tf_regions)
    genome_ragnes = makeGRangesFromDataFrame(tf_regions[,c('chr', 'feature_start','feature_end')])

    matches = as.data.frame( findOverlaps(genome_ragnes, genome_ragnes, minoverlap = 200) )

    matches = matches[matches$queryHits != matches$subjectHits, ]

    str(matches)
     
    if(nrow(matches) > 0){
        #f_one_pair_tf_interaction(match_line, sample_cols, tf_regions)    
        
        overlap_df=data.frame(matches)
        str(overlap_df)
        head(overlap_df)
        overlap_df$name = paste0(tf_regions[overlap_df$queryHits,'feature_tf'], '-', tf_regions[overlap_df$subjectHits,'feature_tf'] )
        overlap_pairs = overlap_df[overlap_df$name %in% valid_interaction$V1,]
        flog.info('%s out of %s is valiad TF interactions',nrow(overlap_pairs), nrow(matches))
        str(overlap_df)
        if(nrow(overlap_pairs)>0){
            tf_interaction_impact = ldply(  apply(overlap_pairs[,1:2], MARGIN = 1, f_one_pair_tf_interaction, sample_cols, tf_regions) )
            dim(tf_interaction_impact)
            head(tf_interaction_impact)
        
            tf_interaction_impact$.id = NULL
            row.names(tf_interaction_impact) = paste0('TF.overlap.',make.names(tf_interaction_impact$feature, unique = TRUE))
            tf_valid_interaction_impact = tf_interaction_impact[tf_interaction_impact$feature %in% valid_interaction$V1,]
            cat('Interaction terms', dim(tf_valid_interaction_impact), '\n')
    
            transcript_data_merge = rbind(transcript_data_tf_concentration, tf_valid_interaction_impact)
            rm(tf_valid_interaction_impact)
            rm(tf_interaction_impact)
        }else{
            cat('Empty overlaps','\n')
        }
        
    }else{
        transcript_data_merge = transcript_data_tf_concentration
    }

    rm(transcript_data_tf_concentration)

    
    ############
    #Promoter - enhancer TF interactions
    ############

    promoter_index=as.numeric(which(tf_regions$type == 'promoter'))
    enhancer_index=as.numeric(which(tf_regions$type == 'enhancer'))
    pair_df=data.frame(promoter = rep(promoter_index, each = length(enhancer_index)), enhancer = rep(enhancer_index, times = length(promoter_index)))
    str(pair_df)
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

    train_data = transcript_data_merge[,sample_cols]
    train_data[is.na(train_data)] = 0

    rm(transcript_data_merge)
    
    train_data_rmdup = train_data[!duplicated(train_data),]
    rm(train_data)
    
    train_data2 = as.data.frame(t(train_data_rmdup))
    
    head10(train_data2)
    if (add_histone == FALSE){
        none_histone_cols = grep('H3K',colnames(train_data2), value = TRUE, invert = TRUE)
        final_train_data = train_data2[, none_histone_cols]
    }else{
        final_train_data  = train_data2
    }

    rm(train_data2)
    rm(train_data_rmdup)
    
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

    
    #final_train_data$population = NULL
    final_train_data[,'population'] = as.character(sample_info[rownames(final_train_data), c('pop')])
    final_train_data[,'gender'] = sample_info[rownames(final_train_data), c('gender')] == 'male'
    final_train_data$population = as.factor(final_train_data$population)
    final_train_data$gender = as.factor(final_train_data$gender)
 
    
    
    additional_cols  =grep('gender|hsa-miR|population', colnames(final_train_data), value = TRUE)
    additional_data =transcript_data[rep('gene.RNASEQ', length(additional_cols)),]
    rownames(additional_data) = additional_cols
    additional_data[,sample_cols] = t(final_train_data[sample_cols, additional_cols])
    
    transcript_data = rbind(transcript_data, additional_data)
    flog.info('After add other feature data')
    #train_model = 'glmnet'
    #train_model = 'enet'
    result <- try(
                        fit  <- f_caret_regression_return_fit(my_train = final_train_data, target_col = 'gene.RNASEQ', learner_name = train_model, tuneGrid),
                        silent=TRUE
                       )
    
    if (class(result) == "try-error"){
        print(result)
        cat('Error in ', transcript_id, '\n')
        next
    }
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
opt_name = f_convet_opts_to_output_dir(opt)
results_dir = f_p('%s/rnaseq/%s/%s/', output_dir, chr_str, opt_name )
dir.create(results_dir)
output_file = f_p('%s/%s.enet',  results_dir, gene)

flog.info('Output file: %s', output_file)

prediction_performance= prediction_performance[-1,]
#prediction_performance$performance = as.numeric(prediction_performance$performance)
#prediction_performance['mean',] =c('mean', mean(prediction_performance[complete.cases(prediction_performance),]$performance))

write.table(prediction_performance, file = output_file , sep = '\t', quote=FALSE, row.names = FALSE)

head(collected_features)
write.table(collected_features, file = f_p('%s.features', output_file), sep = '\t', quote = FALSE, row.names =FALSE)
