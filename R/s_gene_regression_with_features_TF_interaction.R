

setwd('~/expression_var/R/')
source('~/R/s_function.R', chdir = TRUE)
source('s_gene_regression_fun.R')
library(doMC)
library(stringr)
#install.packages('doMC')

library(GenomicRanges)
library("optparse")

option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="dataset file name", metavar="character"),
    make_option(c("--add_histone"), type="character", default='TRUE', help="output file name [default= %default]", metavar="character"),
    make_option(c("--add_miRNA"), type="character", default='TRUE', help="Add miRNA or not", metavar="character"),
    make_option(c("--test"), type="character", default=NULL, help="output file name [default= %default]", metavar="character"),
    make_option(c("--gene"), type="character", default='', help="The name of gene", metavar="character"),
    make_option(c("--model"), type="character", default='enet', help="The machine learning method used, eg. enet, and rfe", metavar="character")
    
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$batch_name)){
    batch_name = '54samples_evalue'
    #batch_name = '462samples_sailfish_quantile'
    add_histone = TRUE
    add_miRNA = TRUE
    test_flag = TRUE
    tuneGrid = NULL
    train_model = 'enet'
    gene=''
}else{
    batch_name = opt$batch_name
    add_histone = opt$add_histone == 'TRUE'
    add_miRNA = opt$add_miRNA == 'TRUE'
    test_flag = opt$test == 'TRUE'
    train_model = opt$model
    gene=opt$model
    cat('Test flag :', test_flag, 'Model ', train_model, '\n')

    tune_list = list( enet = tuneGrid <- expand.grid(.lambda = c(0, 0.001, 0.01, 0.1), .fraction = seq(.05, 1, length = 4 )),
        glm= expand.grid(.alpha = c(0, .1, .2, .4, .6, .8, 1), .lambda = seq(.01, .2, length = 5))
        )

    tuneGrid = tune_list[[train_model]]
    
}


chr_str = 'chr22'
model_str = 'all'
cat('Tune Grid:', is.null(tuneGrid), '\n')

output_dir = f_p('./data/%s/', batch_name)
expression_file = paste0(output_dir, 'rnaseq/gene_regulatory_region_feature_profile.', chr_str )

expression_data = read.table(expression_file, header = TRUE, na.strings = '.')
#head(expression_data)
#tail(expression_data[,1:20])


#cat('Number regulatory regions and empty fragment ids', '\n')
#table(expression_data$type, expression_data$hic_fragment_id == '.')

#table(expression_data$type, expression_data$feature)

#colnames(expression_data)


#Get the expressed transcripts
sample_cols = grep('(NA|HG)[0-9]+', colnames(expression_data),value = T)

cat('The number of samples:', length(sample_cols), '\n')
expression_data = expression_data[,c('chr', 'start', 'end', 'gene','feature','feature_start', 'feature_end', 'hic_fragment_id' ,'type', sample_cols )]
#write.table(sample_cols, './data/output/sample.list', sep = '\t', quote = F, col.names = F)
rna_data = subset(expression_data, feature == 'RNASEQ')
rownames(rna_data) = rna_data$gene
expressed_transcripts = rownames(rna_data)[rowSums(rna_data[, sample_cols]) > 0]
cat('The number of expressed transcripts:', length(expressed_transcripts), '\n')

#head(rna_data)
#str(rna_data)
#Get the genes with at least one elements
gene_with_regulatory_elements = as.data.frame.matrix(table(expression_data$gene, expression_data$feature))
head(gene_with_regulatory_elements)
regulated_genes_names = rownames(gene_with_regulatory_elements)[ rowSums(gene_with_regulatory_elements) > 10]
cat('The number of regulated transcripts:', length(regulated_genes_names), '\n')

#Get the inverstigeed genes
genes_names = intersect(expressed_transcripts, regulated_genes_names)
cat('The number of investigated transcripts:', length(genes_names), '\n')
0

length(genes_names)

#print(genes_names[1:20])

prediction_performance = data.frame(gene = 'mean', performance = '0', SD = '0', stringsAsFactors = F)
collected_features = data.frame()
sample_info = read.table(f_p('%s/chr_vcf_files/integrated_call_samples_v3.20130502.ALL.panel', output_dir ), header = TRUE, row.names = 1)




#Read the miRNA interaction and expression data.
miRNA_target_table = read.table('./data/raw_data/miRNA/miRNA_ensemble.txt', header = TRUE)
miRNA_expression = read.table('./data/raw_data/miRNA/GD452.MirnaQuantCount.1.2N.50FN.samplename.resk10.txt', header = TRUE)

i = 10

sample_cols = intersect(sample_cols, colnames(miRNA_expression))
#length(sample_cols)
non_sample_cols = setdiff(colnames(expression_data), sample_cols)
print(non_sample_cols)

#Read the TF expression data.
tf_gene_id = read.table('./data/raw_data/rnaseq/tf_ensemble_id.txt', header = T, stringsAsFactors = FALSE)
#tf_gene_id
rownames(tf_gene_id) = tf_gene_id$tf
expression_data$feature_tf =''
str(expression_data)
#head(expression_data)

tf_gene_id['P300',]

expression_data$feature_tf = (tf_gene_id[expression_data$feature,'external_name'])


#head(tf_gene_id)
gene_expression = read.table(f_p('%s/rnaseq/transcript_data.bed', output_dir), header = T, stringsAsFactors = FALSE)
rownames(gene_expression) = gene_expression$transcript_id
#str(tf_gene_id)

tf_gene_expression = gene_expression[tf_gene_id$ensembl_gene_id,]
#dim(tf_gene_expression)
#head(tf_gene_expression[,1:20])

rownames(tf_gene_expression) = tf_gene_id$tf


#head(tf_gene_expression,1)

range(tf_gene_expression[1, sample_cols])

valid_interaction = read.table('./data/raw_data/biogrid/tf_interactions.txt')

#head(valid_interaction)

sum(duplicated(expression_data[, non_sample_cols[1:8]]))

for (i in 1:length(genes_names)){
    if(i == 3 & test_flag == TRUE){
        print( 'First 500 genes')
        break
    }
    
    transcript_id = genes_names[i]
    cat('\n','============',i, transcript_id, '==============','\n')
    sum(gene_with_regulatory_elements[transcript_id,])
    gene_with_regulatory_elements[transcript_id,]
    transcript_data = expression_data[expression_data$gene == transcript_id,]
    dim(transcript_data)
    head(transcript_data[,1:15])
    transcript_data = transcript_data[!duplicated(transcript_data),]
    head(transcript_data)
    table(transcript_data$feature)
    rownames(transcript_data) = make.names(paste0(transcript_data$type, '-',transcript_data$feature), unique = T)
    transcript_data[is.na(transcript_data)]=0 #Some TF regions don't have variations 
    
    #######Add the TF concentration data#########
    tmp=tf_gene_expression[as.character(transcript_data$feature),sample_cols]
    tmp[is.na(tmp)]=1 #Set the non-TF rows to 1
    scaled_tmp=t(apply(tmp, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))) + 1)
    scaled_tmp[is.na(scaled_tmp)] = 1 #For the DNASE regions.

    head(tf_gene_expression)
    head(scaled_tmp[,1:15], 1000)
    dim(tmp)
    
    transcript_data = transcript_data[!duplicated(transcript_data[, non_sample_cols[1:8]]),]
    dim(transcript_data)
    
    subset(transcript_data[, non_sample_cols], feature == 'TBP')
    subset(transcript_data[, non_sample_cols], hic_fragment_id == '24949953')
    transcript_data_tf_concentration = transcript_data
    transcript_data_tf_concentration[, sample_cols] = transcript_data[, sample_cols] * scaled_tmp[, sample_cols]
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

    if(nrow(matches) > 0){
        #f_one_pair_tf_interaction(match_line, sample_cols, tf_regions)    
        library(plyr)
        tf_interaction_impact = ldply(  apply(matches, MARGIN = 1, f_one_pair_tf_interaction, sample_cols, tf_regions) )
        dim(tf_interaction_impact)
        head(tf_interaction_impact)


        tf_interaction_impact$.id = NULL
        row.names(tf_interaction_impact) = paste0('TF.overlap.',make.names(tf_interaction_impact$feature, unique = TRUE))
        tf_valid_interaction_impact = tf_interaction_impact[tf_interaction_impact$feature %in% valid_interaction$V1,]
        cat('Interaction terms', dim(tf_valid_interaction_impact), '\n')
    
        transcript_data_merge = rbind(transcript_data_tf_concentration, tf_valid_interaction_impact)
    }else{
        transcript_data_merge = transcript_data_tf_concentration
    }

    

    
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
        #str(promoter_pairs)

        if(nrow(promoter_pairs) > 0){
            
            promoter_interaction_impact = ldply(  apply(promoter_pairs[,1:2], MARGIN = 1, f_one_pair_tf_interaction, sample_cols, tf_regions) )
            dim(promoter_interaction_impact)
            promoter_interaction_impact$.id=NULL
            rownames(promoter_interaction_impact) = make.names(paste0('P_E.', promoter_interaction_impact$feature), unique = TRUE)
            
            dim(transcript_data_merge)
            transcript_data_merge = rbind(transcript_data_merge, promoter_interaction_impact)
        }else{
            cat('Empty promoter-enhancers!', '\n')
        }
    }


    rownames(transcript_data_merge)
    train_data = transcript_data_merge[,sample_cols]
    train_data[is.na(train_data)] = 0
    
    train_data_rmdup = train_data[!duplicated(train_data),]

    train_data2 = as.data.frame(t(train_data_rmdup))
    #head(train_data2)
    if (add_histone == FALSE){
        none_histone_cols = grep('H3K',colnames(train_data2), value = TRUE, invert = TRUE)
        final_train_data = train_data2[, none_histone_cols]
    }else{
        final_train_data  = train_data2
    }
 
    
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
    train_model = 'enet'
    
    result <- try(
                        fit  <- f_caret_regression_return_fit(my_train = final_train_data, target_col = 'gene.RNASEQ', learner_name = train_model, tuneGrid),
                        silent=TRUE
                       )
    if (class(result) == "try-error"){
        print(result)
        cat('Error in ', transcript_id, '\n')
        next
    }


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
    
    #fit$bestTune
    
    #fit  <- f_caret_regression_return_fit(my_train = final_train_data, target_col = 'gene.RNASEQ', learner_name = 'glmnet')
    
    #head(final_train_data)
    #print(fit$results)
    max_performance = fit$max_performance
    RsquaredSD =fit$results[rownames(fit$bestTune), 'RsquaredSD' ]
    
    cat("\n",'performance', max_performance, 'SD', RsquaredSD, 'Number of features:', nrow(features_df) )
    prediction_performance = rbind(prediction_performance, c(as.character(transcript_id), as.character(max_performance), as.character(RsquaredSD)))
    #print(prediction_performance)
}

output_file = f_p('%s/output/lasso.prediction.%s.%s.%s.%s.%s.txt', output_dir, chr_str, model_str, ifelse(add_histone,  'addHis', 'rmvHis'), ifelse(add_miRNA, 'addMiR', 'rmvMiR'), train_model )

prediction_performance= prediction_performance[-1,]
prediction_performance$performance = as.numeric(prediction_performance$performance)
prediction_performance['mean',] =c('mean', mean(prediction_performance[complete.cases(prediction_performance),]$performance), '0')

write.table(prediction_performance, file = output_file , sep = '\t', quote=FALSE, row.names = TRUE)

head(collected_features)
write.table(collected_features, file = f_p('%s.features', output_file), sep = '\t', quote = FALSE, row.names =FALSE)

performance_data = read.table(file = output_file, sep = '\t')
cat("Mean R-square:",  mean(performance_data[complete.cases(performance_data),]$performance), '\n')









