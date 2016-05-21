
#This script combine the results from s_regression_for_one_gene.R

#Sync the data from clustdell.

library(matrixStats)
setwd('~/expression_var/R/')
source('~/R/s_function.R', chdir = TRUE)
source('s_gene_regression_fun.R')
source('s_summary_fun.R')
library(stringr)
#install.packages('doMC')
library(dplyr)
library("optparse")

option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="dataset file name", metavar="character"),
    make_option(c("--chr"),      type="character", default='chr22', help="chromosome name", metavar="character")    
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$batch_name)){
    #batch_name = '54samples_evalue'
    #batch_name = '462samples_sailfish_quantile'
    #batch_name = '462samples_snyder_norm'
    #batch_name = '462samples_quantile'
    #batch_name = '462samples_quantile_rmNA'
    batch_name = '462samples_log_quantile'
    #chr_str = 'chr22'
    chr_str = 'chr22'
    
}else{
    batch_name = opt$batch_name
    chr_str = opt$chr
}

##Compare remove and add gm12878 bias.

mode_list = c('rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_rm.gm12878',
              'rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878')

small_sample_performance = f_collect_performance_in_mode_list(batch_name, chr_str, mode_list, rsync_flag = TRUE)
colnames(small_sample_performance) = c('rm.gm12878', 'add.gm12878')
head(small_sample_performance)
apply(small_sample_performance, 2, mean)


###Compare the small samples.
###Conclusion: enet is better than gmlnet in small samples.
###29 samples gives better accuracy than 54 smaples.
mode_list = c('add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.TF.exp.only_rm.predict.TF_add.YRI_population.29YRI_TF.exp.type.fakeTF_add.gm12878.TURE',
              'add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.TF.exp.only_rm.predict.TF_add.YRI_population.54snyder_TF.exp.type.fakeTF_add.gm12878.TURE',
              'add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.TF.exp.only_rm.predict.TF_add.YRI_population.29snyder_TF.exp.type.fakeTF_add.gm12878.TURE',
              'rm.histone_rm.miRNA_rm.TF.exp_model.enet_rm.TF.exp.only_rm.predict.TF_add.YRI_population.54snyder_TF.exp.type.fakeTF_add.gm12878')

small_sample_performance = f_collect_performance_in_mode_list(batch_name, chr_str, mode_list, rsync_flag = TRUE)
colnames(small_sample_performance) = c('glmnet.29YRI', 'glmnet.54snyder', 'glmnet.29snyder', 'enet.54snyder')
head(small_sample_performance)
apply(small_sample_performance, 2, mean)


stop()

###########Compare the 54samples enet and glmnet##########################
loc_batch_name = '462samples_log_quantile'
mode_list = c('rm.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.TF.exp.only_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878',
              'rm.histone_rm.miRNA_rm.TF.exp_model.enet_rm.TF.exp.only_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878')
sample54_performance = f_collect_performance_in_mode_list(loc_batch_name, chr_str, mode_list, rsync_flag = TRUE)

head(sample54_performance)
colnames(sample54_performance) = c('glmnet', 'enet')

f_plot_two_columns(sample54_performance, c('glmnet', 'enet'))



for(loc_mode in mode_list){
    performance = f_summary_regression_results( loc_batch_name, chr_str, loc_mode, rsync_flag = FALSE)
    mean_stats <- performance %>% filter(performance > 0.25) %>% summarise(performance = mean(performance), SD = mean(SD), selected_feature = mean(seleted_features))
    print(mean_stats)
}


shared_rows = intersect(rownames(sample54_performance), rownames(small_sample_performance))
length(shared_rows)
combined_data = cbind(small_sample_performance[shared_rows,2], sample54_performance[shared_rows,])


head10(collected_performance)

mode_list =c('add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_rm.predict.TF_rm.YRI')
collected_performance = f_collect_performance_in_mode_list(batch_name, chr_str, mode_list, rsync_flag = TRUE)

stop()
###########Compare the random gene expression 462samples_random###########################

mode_list = c('add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.TF.exp.only_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF')
collected_performance = f_collect_performance_in_mode_list('462samples_random', chr_str, mode_list, rsync_flag = TRUE)

random_performance = f_summary_regression_results('462samples_random', chr_str, mode_list[1])
apply(random_performance[,c('performance', 'SD')], 2, mean)

norm_performance = f_summary_regression_results('462samples_log_quantile', chr_str, 'add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_rm.predict.TF_add.YRI')

random_performance %>% filter(performance > 0.25) %>% summarise(performance = mean(performance), SD = mean(SD))
norm_performance %>% filter(performance > 0.25) %>% summarise(performance = mean(performance), SD = mean(SD))





stop()

#########Compare add miRNA or not#####################
mode_list =c('add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_rm.predict.TF_add.YRI',
    'add.histone_add.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_rm.predict.TF_add.YRI_population.all')
                     

collected_performance = f_collect_performance_in_mode_list(batch_name, chr_str, mode_list, rsync_flag = TRUE)
colnames(collected_performance) = c('rm_miRNA', 'add_miRNA')
collected_performance = as.data.frame(collected_performance)
ggplot((collected_performance), aes(y = add_miRNA, x = rm_miRNA )) + geom_point() + geom_abline()



stop()

#####Compare three machine learning methods#################

mode_list =c('add.histone_rm.miRNA_rm.TF.exp_model.rf_rm.permutation_rm.TF.exp.only_rm.predict.TF_add.YRI_population.all',
                     'add.histone_rm.miRNA_rm.TF.exp_model.gbm_rm.permutation_rm.TF.exp.only_rm.predict.TF_add.YRI_population.all',
                     'add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_rm.predict.TF_add.YRI')

collected_performance = f_collect_performance_in_mode_list(batch_name, chr_str, mode_list, rsync_flag = FALSE)
colnames(collected_performance) = c('sequence_prediction', 'add_TF_expression')
collected_performance = as.data.frame(collected_performance)
ggplot((collected_performance), aes(y = add_TF_expression, x = sequence_prediction )) + geom_point() + geom_abline()


stats_collected = data.frame()
#Read gene expression
rnaseq_data=read.table('./data/462samples_sailfish/rnaseq/GEUVADIS.Gene.DATA_MATRIX', header = T, sep = ' ')
rownames(rnaseq_data) = rnaseq_data$gene
rnaseq_data$gene = NULL
head10(rnaseq_data)
for (loc_mode in mode_list){
    performance_df = f_summary_regression_results(batch_name, chr_str, loc_mode, rsync_flag = FALSE)
    accurate_genes = rownames(subset(performance_df, performance > 0.25))
    low_genes = rownames(subset(performance_df, performance <= 0.25))
    stats_one_model <- performance_df %>% filter(performance > 0.25) %>% dplyr::summarise(number = length(performance), performance = mean(performance_df$performance) ,num_feature = mean(num_feature), selected_features = mean(seleted_features))
    stats_one_model[, c('high_mean', 'high_std', 'low_mean', 'low_std')] = c(mean(rowMeans(rnaseq_data[accurate_genes,])),
    mean(rowSds(as.matrix(rnaseq_data[low_genes,]))),
    mean(rowSds(as.matrix(rnaseq_data[accurate_genes,]))),
    mean(rowMeans(rnaseq_data[low_genes,]))
    )
    stats_collected = rbind(stats_collected, stats_one_model)
}

hist(performance_df$performance, main='Performance of Elastic Net', xlab = 'R-square')

print(stats_collected)
head10(performance_df)


###########Feature stats#############
loc_mode = 'add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_rm.predict.TF_add.YRI'
features_df = f_summary_regression_results(batch_name, chr_str, loc_mode, rsync_flag = FALSE, return_features = TRUE)
performance_df = f_summary_regression_results(batch_name, chr_str, loc_mode, rsync_flag = FALSE, return_features = FALSE)
accurate_genes = rownames(subset(performance_df, performance > 0.25))


library(plyr)
features_df[, c('gene', 'feature')] = ldply(str_split(features_df$name,'[|]'))

accurate_features = features_df[features_df$gene %in% accurate_genes,]
dim(accurate_features)

write.table(sort(table(str_replace(accurate_features$feature, '[.][0-9]*$', '')), decreasing = TRUE), 'feature.table', quote = FALSE, sep = '\t')






stop()

####Compare the add YRI and removing YRI.

mode_list =c('add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_rm.predict.TF_rm.YRI',
                     'add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_rm.predict.TF_add.YRI')

collected_performance = f_collect_performance_in_mode_list(batch_name, chr_str, mode_list, rsync_flag = TRUE)
colnames(collected_performance) = c('sequence_prediction', 'add_TF_expression')
collected_performance = as.data.frame(collected_performance)
ggplot((collected_performance), aes(y = add_TF_expression, x = sequence_prediction )) + geom_point() + geom_abline()




stop()
###########################

batch_list = c('462samples_quantile', '462samples_quantile_rmNA', '462samples_var_quantile', '462samples_snyder_norm')
for (loc_batch_name in batch_list){
    cat('======', batch_name, '========\n')
    #performance = f_summary_regression_results(loc_batch_name, chr_str, mode_name)
    Sys.sleep(1)
}


mode_list = c('add.histone_add.miRNA_add.TF.exp_model.glmnet_add.permutation',
    'add.histone_add.miRNA_add.TF.exp_model.glmnet_rm.permutation',
    'add.histone_add.miRNA_rm.TF.exp_model.glmnet_rm.permutation',
      'add.histone_rm.miRNA_add.TF.exp_model.glmnet_rm.permutation',
    'add.histone_add.miRNA_add.TF.exp_model.glmnet_add.permutation_add.TF.exp.only'
    )

#Compare GBM and glmnet
mode_list =c('add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_rm.predict.TF',
                     'add.histone_rm.miRNA_rm.TF.exp_model.gbm_rm.permutation_rm.TF.exp.only_rm.predict.TF')

collected_performance = f_collect_performance_in_mode_list(batch_name, chr_str, mode_list, rsync_flag = TRUE)
colnames(collected_performance) = c('glmnet', 'gbm')
collected_performance = as.data.frame(collected_performance)
ggplot((collected_performance), aes(y = glmnet, x = gbm )) + geom_point() + geom_abline()


#Compare the difference
mode_list =c('add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_rm.predict.TF',
                    'add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_add.TF.exp.only_rm.predict.TF')

collected_performance = f_collect_performance_in_mode_list(batch_name, chr_str, mode_list, rsync_flag = TRUE)
colnames(collected_performance) = c('sequence_prediction', 'add_TF_expression')
collected_performance = as.data.frame(collected_performance)
ggplot((collected_performance), aes(y = add_TF_expression, x = sequence_prediction )) + geom_point() + geom_abline()



#Compare the raw and raw + TF prediction
mode_list = c('add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_add.predict.TF',
                      'add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only')
collected_performance = f_collect_performance_in_mode_list(batch_name, chr_str, mode_list, rsync_flag = FALSE)
colnames(collected_performance) = c('add_TF_prediction', 'sequence_prediction')
collected_performance = as.data.frame(collected_performance)
ggplot((collected_performance), aes(y = add_TF_prediction, x = sequence_prediction )) + geom_point() + geom_abline()


head10(collected_performance)
collected_performance = f_sort_by_col(collected_performance, 'add_TF_prediction', decr_flag = TRUE)
rownames(collected_performance[1:20,])

write.table(collected_performance[1:20,], f_p('./data/%s/rnaseq/%s/top20', batch_name, chr_str), quote=FALSE)


stop()
##########################
###Get the bed files for high accuracy and low accuracy genes.
#######################

performance_threshold = 0.1

head(performance_df)

all.entrezgene = read.table('./data/raw_data/rnaseq/all.ensemble.genes',sep='\t', header = TRUE)
row.names(all.entrezgene) = all.entrezgene$ensembl_transcript_id

head(all.entrezgene)
duplicated_rows = duplicated(all.entrezgene$ensembl_gene_id)
all.entrezgene = all.entrezgene[!duplicated_rows,]
row.names(all.entrezgene) = all.entrezgene$ensembl_gene_id

#Change to the promoter regions.
head(all.entrezgene)
rna_seq = all.entrezgene[, c('chromosome_name','transcript_start', 'transcript_end', 'strand')]
rna_seq = rna_seq[complete.cases(rna_seq),]
head(rna_seq)
colnames(rna_seq) = c('chr','start','end', 'strand')
tmp_data = rna_seq[,c('chr','start','end', 'strand')]
forward_strand = rna_seq$strand == 1
reverse_strand  = rna_seq$strand == -1
cat('Forward strand:', sum(forward_strand), 'Reverse ', sum(!forward_strand) ,'\n')

rna_seq[forward_strand,'start'] = tmp_data[forward_strand, 'start'] - 2000
rna_seq[forward_strand,'end'] = tmp_data[forward_strand, 'start'] + 2000
rna_seq[reverse_strand,'start'] = tmp_data[reverse_strand, 'end'] - 2000
rna_seq[reverse_strand,'end'] = tmp_data[reverse_strand, 'end'] + 2000

library(plyr)
gene_ids=ldply(str_split(performance_df$gene, '[.]'))
performance_df[,c('chr','start','end', 'strand')] = rna_seq[gene_ids$V1,]
performance_df = performance_df[complete.cases(performance_df),]
head(performance_df)

performance_threshold = 0.1
table(performance_df$performance > performance_threshold)
performance_df$chr = paste0('chr', performance_df$chr)

high_group = subset(performance_df, performance >= performance_threshold )
low_group = subset(performance_df, performance < performance_threshold )


write.table(high_group[, c('chr','start','end', 'strand')], f_p('%s/high.bed', output_dir), quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
write.table( performance_df[, c('chr','start','end', 'strand')], f_p('%s/low.bed', output_dir), quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)










