
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
source('~/R/s_ggplot2_theme.R', chdir = TRUE)


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



#loc_mode = mode_list[1]
f_summary_selected_features <- function(loc_batch_name, chr_list, loc_mode, output_feature_file, performance_threshold){
    
    accurate_features = data.frame()
    for (chr_str in chr_list){
        features_df = f_summary_regression_results(loc_batch_name, chr_str, loc_mode, rsync_flag = FALSE, return_features = TRUE)
        performance_df = f_summary_regression_results(loc_batch_name, chr_str, loc_mode, rsync_flag = FALSE, return_features = FALSE)
        accurate_genes = rownames(subset(performance_df, performance > performance_threshold))

        library(plyr)
        features_df[, c('gene', 'feature')] = ldply(str_split(features_df$name,'[|]'))

        accurate_features = rbind( accurate_features, features_df[features_df$gene %in% accurate_genes,])
    }
    dim(accurate_features)
    feature_rank_table = sort(table(str_replace(accurate_features$feature, '[.][0-9]*$|[.]rs[0-9]*', '')), decreasing = TRUE)
    write.table(feature_rank_table, output_feature_file, quote = FALSE, sep = '\t')
    return (feature_rank_table)
}

mode_list = c('rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_add.batch.random',
              'rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_rm.batch.random')

mode_list =c('rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.None_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_add.batch.random',
             'rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.None_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_rm.batch.random')

non_pop_list = c(All='rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.all',
                 SNP='rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.SNP',
                 TF='rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.TF',
                 random = 'rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.random')


pop_list = c(All = 'rm.histone_rm.miRNA_model.glmnet_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.all',
             SNP = 'rm.histone_rm.miRNA_model.glmnet_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.SNP',
             random = 'rm.histone_rm.miRNA_model.glmnet_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.random',
             TF = 'rm.histone_rm.miRNA_model.glmnet_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.TF')


genebody_list = c(SNP='rm.histone_rm.miRNA_model.glmnet_population.all_add.gm12878_new.batch.462samples.genebody_batch.mode.SNP',
                            TF = 'rm.histone_rm.miRNA_model.glmnet_population.all_add.gm12878_new.batch.462samples.genebody_batch.mode.TF')


####################Compare the 462samples peer#####################
#chr_num_list = c(22, 2, 10, 20)
chr_num_list = c(22)


sample_performance_merge = data.frame()
for (chr_num in chr_num_list){
    loc_batch_name = '462samples_quantile_rmNA'
    chr_str = paste0('chr', chr_num)
    dir.create(f_p('./data/%s/rnaseq/%s/', loc_batch_name, chr_str))
    

    mode_list = genebody_list
    sample_performance = f_collect_performance_in_mode_list(loc_batch_name, chr_str, mode_list, rsync_flag = TRUE)
    sample_performance = as.data.frame(sample_performance)
    
    colnames(sample_performance)
    head(sample_performance)
    sample_performance$chr = chr_str
    sample_performance_merge = rbind(sample_performance_merge, sample_performance)
    
}
results_dir = f_p('./data/%s/rnaseq/results/', loc_batch_name)

head(sample_performance_merge)
str(sample_performance_merge)
table(sample_performance_merge$mode)

library(reshape2)
library(futile.logger)
library(dplyr)
sample_performance_merge %>% group_by(mode) %>%
    dplyr::summarise(good_predictions = sum(performance > 0.1),
                               mean_performance = mean(performance),
                               mean_SD = mean(SD),
                               feature_num = mean(selected_features),
                     max_feature = max(selected_features))

table(sample_performance_merge$mode)

threshold = quantile( subset(sample_performance_merge, mode == 'random')$performance , probs = 0.96 )



good_predictions <- unique(as.character((sample_performance_merge %>% filter(performance > threshold, performance < 1 ,mode != 'random') )$gene))

sample_performance_merge %>% filter(selected_features > 400)

basic_stats <- sample_performance_merge %>% filter(performance > threshold) %>% group_by(mode) %>%
    dplyr::summarise(Total = length(good_predictions),
                               good_predictions = sum(performance > threshold),
                               mean_performance = mean(performance),
                               #mean_SD = mean(SD),
                               mean_input_feature = mean(num_feature), 
                               mean_selected_features = mean(selected_features),
                               max_feature = max(selected_features))

basic_stats <- sample_performance_merge %>% filter(gene %in% good_predictions) %>% group_by(mode) %>%
    dplyr::summarise(Total = length(good_predictions),
                               good_predictions = sum(performance > threshold),
                               mean_performance = mean(performance),
                               #mean_SD = mean(SD),
                               mean_input_feature = mean(num_feature), 
                               mean_selected_features = mean(selected_features),
                               max_feature = max(selected_features))

write.table(format(as.data.frame(basic_stats), digits = 3), f_p('%s/basic_stats.txt', results_dir), quote = FALSE, sep = '\t', row.names = FALSE)

head(sample_performance_merge)
library(limma)


####Venn plot######


TF <- subset( sample_performance_merge, gene  %in% good_predictions & mode == 'TF')$performance >= threshold
SNP<- subset( sample_performance_merge, gene %in% good_predictions & mode == 'SNP')$performance >= threshold
All<-subset( sample_performance_merge, gene %in% good_predictions & mode == 'All')$performance >= threshold
c3 <- cbind(TF, SNP, All)
a <- vennCounts(c3)
vennDiagram(a)


good_performance <- sample_performance_merge %>% filter(gene %in% good_predictions)
qplot(subset(good_performance, mode == 'SNP')$performance, subset(good_performance, mode == 'All')$performance) + geom_abline()
wilcox.test(subset(good_performance, mode == 'SNP')$performance, subset(good_performance, mode == 'All')$performance, paired = T)
qplot(subset(good_performance, mode == 'All')$selected_features, subset(good_performance, mode == 'SNP')$selected_features) +
    geom_abline() + coord_fixed()

wilcox.test(subset(good_performance, mode == 'SNP')$selected_features, subset(good_performance, mode == 'TF')$selected_features, paired = T)
qplot(subset(good_performance, mode == 'TF')$performance, subset(good_performance, mode == 'All')$performance) +
    geom_abline() + coord_fixed()
    



length(unique(sample_performance$gene))
sample_performance_merge$Model = str_replace(sample_performance_merge$mode, 'random', 'Random')


four_modes_plot <-ggplot(sample_performance_merge, aes(x = performance , color = mode )) + geom_density()
TF_random_plot <- ggplot(subset(sample_performance_merge, mode != 'All' & mode != 'SNP'),
                         aes(x = performance, color = mode )) + geom_density()
ggsave(f_p('%s/TF_random_performance.tiff', results_dir), plot = TF_random_plot)


TF_SNP_All_performance_plot <- ggplot(subset(sample_performance_merge, mode != 'random'),
                                      aes(x = performance, color = mode )) + geom_density() + theme_set(theme_grey(base_size = 22)) 
ggsave(f_p('%s/TF_SNP_All_performance.tiff', results_dir), plot = TF_SNP_All_performance_plot)


All_performance_plot <- ggplot(subset(sample_performance_merge, performance >0),
                                      aes(x = performance, color = Model )) + xlab('Adjusted R square') +
    geom_density() + theme_set(theme_grey(base_size = 18))
                                                                                                         
All_performance_plot
ggsave(f_p('%s/All_performance.tiff', results_dir), plot = All_performance_plot)



stop()
#####Get high accracte data#####



#####Count the feature categories######
feature_table = f_summary_selected_features(loc_batch_name, c('chr2', 'chr20', 'chr10', 'chr22'), mode_list[1], output_feature_file = 'feature.table', threshold )
feature_table[1]/sum(feature_table)
head(feature_table)


tf_feature_table = f_summary_selected_features(loc_batch_name, c('chr2', 'chr20', 'chr10', 'chr22'), mode_list[3], output_feature_file = 'feature.table.tf', threshold )



#######Find top features####################
all_performance <- sample_performance_merge %>% filter(mode == 'All', chr == 'chr22')

(f_sort_by_col(all_performance, 'performance', decr_flag = T)[1:5,])$gene




##############Adjust R-square#####################
good_features = subset(small_features, adjust.rsq>0.05 )
table(good_features$mode)

adjust_basic_stats <- good_features %>% group_by(mode) %>%
    dplyr::summarise(good_predictions = length(mode),
                               mean_performance = mean(adjust.rsq),
                               feature_num = mean(selected_features),
                                max_feature = max(selected_features))



##############Sort features#################

target_gene = 'ENSG00000172404.4'
feature_file = f_p('./data/462samples_quantile_rmNA/rnaseq/chr22/%s/%s.enet.features', mode_list[1], target_gene)
feature_table = read.table(feature_file, header = T)
feature_table$abs_score = abs(feature_table$score)
sort_table = f_sort_by_col(feature_table, index = 'abs_score', T)
sort_table$abs_score = NULL
head(sort_table)
sort_table$chr ='chr22'
write.table(sort_table[,c('name', 'score')], file = f_p('%s.sort', feature_file), quote = F, sep = '\t', row.names = F)

feature_file = f_p('./data/462samples_quantile_rmNA/rnaseq/chr22/%s/%s.enet.features', mode_list[3], target_gene)
feature_table = read.table(feature_file, header = T)
feature_table$abs_score = abs(feature_table$score)
sort_table = f_sort_by_col(feature_table, index = 'abs_score', T)
sort_table$abs_score = NULL
head(sort_table)
sort_table$chr ='chr22'
write.table(sort_table[,c('name', 'score')], file = f_p('%s.sort', feature_file), quote = F, sep = '\t', row.names = F)
