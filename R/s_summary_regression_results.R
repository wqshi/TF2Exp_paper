
#This script combine the results from s_regression_for_one_gene.R
#Sync the data from clustdell.
library(reshape2)
library(futile.logger)
library(dplyr)
library(matrixStats)
setwd('~/expression_var/R/')
source('~/R/s_function.R', chdir = TRUE)
source('s_gene_regression_fun.R')
source('s_summary_fun.R')
library(stringr)
library("optparse")
source('~/R/s_ggplot2_theme.R', chdir = TRUE)
options(dplyr.width = Inf)

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
mode_list = c('rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_add.batch.random',
              'rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_rm.batch.random')

mode_list =c('rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.None_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_add.batch.random',
             'rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.None_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_rm.batch.random')

non_pop_list = c(All='rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.all',
                 SNP='rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.SNP',
                 TF='rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.TF',
                 random = 'rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.random')


pop_list = c(
    AlltfShuffle ='rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.AlltfShuffle',
    noInteract = 'rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.noInteract',
    All             = 'rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.all',
             SNP = 'rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.SNP',
    fakeInteract = 'rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.fakeInteract',
             #random = 'rm.histone_rm.miRNA_model.glmnet_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.random',
             TF = 'rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.TF')

genebody_list = c( SNP='rm.histone_rm.miRNA_model.glmnet_population.all_add.gm12878_new.batch.462samples.genebody_batch.mode.SNP',
                            TF = 'rm.histone_rm.miRNA_model.glmnet_population.all_add.gm12878_new.batch.462samples.genebody_batch.mode.TF',
                            SNPinTF = 'rm.histone_rm.miRNA_model.glmnet_population.all_add.gm12878_new.batch.462samples.genebody_batch.mode.SNPinTF')

cvglmnet_list = c(
    AlltfShuffle ='rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.AlltfShuffle',
    #noInteract = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.noInteract',
    All             = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.all',
    #AllnoInteract = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.AllnoInteract',
    SNP = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.SNP'
    #fakeInteract = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.fakeInteract',
             #random = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.random',
             #TF = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.TF'
)

penalty_list = c(TF = 'rm.histone_rm.miRNA_model.glmnet_add.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.TF',
                 all = 'rm.histone_rm.miRNA_model.glmnet_add.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.all')

loc_batch_name = '462samples_quantile_rmNA'
mode_list = cvglmnet_list


#loc_batch_name = '462samples_genebody'
#mode_list = genebody_list
chr_num_list = c(22)
####################Compare the 462samples peer#####################
#chr_num_list = c(22, 2, 10, 20)

sample_performance_merge = f_collect_performance_in_multiple_chrs(loc_batch_name, mode_list, chr_num_list)
sample_performance_merge$performance = sample_performance_merge$rsq
results_dir = f_p('./data/%s/rnaseq/results/', loc_batch_name)


table(sample_performance_merge$mode)

#ggplot(sample_performance_merge, aes(x=train_performance, y = performance, color = mode)) + geom_point() + geom_abline()
#ggplot(sample_performance_merge, aes(x=train_performance - performance, color = mode)) + geom_density()


sample_performance_merge %>% group_by(mode) %>%
    dplyr::summarise(good_predictions = sum(performance > 0.1),
                               mean_performance = mean(performance),
                               mean_SD = mean(SD),
                               feature_num = mean(selected_features),
                     max_feature = max(selected_features))
head(sample_performance_merge)
sample_performance_merge %>% filter(gene == 'ENSG00000239900.7')


sample_performance_merge %>% filter( performance > 0.5 * train_performance ) %>% group_by(mode) %>%
    dplyr::summarise(good_predictions = sum(performance > 0.1),
                               mean_performance = mean(performance),
                               mean_SD = mean(SD),
                               feature_num = mean(selected_features),
                     max_feature = max(selected_features))


gene_count = table(sample_performance_merge$gene) == length(mode_list)
intersect_genes = names(gene_count)[gene_count]
union_good_genes <- unique(as.character((sample_performance_merge %>% filter(performance > threshold, performance < 1 ,mode != 'random') )$gene))



#threshold = quantile( subset(sample_performance_merge, mode == 'random')$performance , probs = 0.99 )
threshold = 0.05



intersect_stats <- sample_performance_merge %>% filter(performance > threshold, gene %in% intersect_genes) %>% group_by(mode) %>%
    dplyr::summarise(Total = length(intersect_genes),
                               good_predictions = sum(performance > threshold),
                               mean_performance = mean(performance),
                               #mean_SD = mean(SD),
                               mean_input_feature = mean(num_feature), 
                               mean_selected_features = mean(selected_features),
                               max_feature = max(selected_features))
print(intersect_stats)
write.table(format(as.data.frame(basic_stats), digits = 3), f_p('%s/basic_stats.txt', results_dir), quote = FALSE, sep = '\t', row.names = FALSE)



####Venn plot######
library(limma)

TF <- subset( sample_performance_merge, gene  %in% intersect_genes & mode == 'TF')$performance >= threshold
SNP<- subset( sample_performance_merge, gene %in% intersect_genes & mode == 'SNP')$performance >= threshold
All<-subset( sample_performance_merge, gene %in% intersect_genes & mode == 'All')$performance >= threshold
c3 <- cbind(TF, SNP, All)
a <- vennCounts(c3)
vennDiagram(a)

table(sample_performance_merge$mode)

good_performance <- sample_performance_merge %>% filter(gene %in% intersect_genes, rsq > 0.05)
dim(sample_performance_merge)
dim(good_performance)



good_performance2 = good_performance %>% filter(mode == 'SNP' | mode == 'All', performance > 0.05)
dim(good_performance2)
table(good_performance2$mode)
ggplot(good_performance2, aes(x = performance , color = mode )) + geom_density()
ggplot(, aes(x = performance , color = mode )) + geom_density()



f_plot_performance_and_stats_test(good_performance, 'All', 'SNP' )
f_plot_performance_and_stats_test(good_performance, 'SNP', 'AllnoInteract' )
f_plot_performance_and_stats_test(good_performance, 'SNP', 'AlltfShuffle')
f_plot_performance_and_stats_test(good_performance, 'All', 'AlltfShuffle' )
#f_plot_performance_and_stats_test(good_performance, 'All', 'TF')
#f_plot_performance_and_stats_test(good_performance, 'TF', 'noInteract')
#f_plot_performance_and_stats_test(good_performance, 'TF', 'fakeInteract')
#f_plot_performance_and_stats_test(good_performance, 'noInteract', 'fakeInteract') #Not significant

#good_intersect genes.
sample_performance_merge$gene = as.character(sample_performance_merge$gene)
robust_performance <- sample_performance_merge %>% filter( rsq > 0.5 * train_performance, performance > 0.05) 
robust_genes = names(gene_count)[gene_count]
dim(robust_performance)
f_plot_performance_and_stats_test(robust_performance, 'SNP', 'All')
f_plot_performance_and_stats_test(robust_performance, 'TF', 'All')
f_plot_performance_and_stats_test(robust_performance, 'TF', 'SNP')
f_plot_performance_and_stats_test(robust_performance, 'TF', 'noInteract')

robust_performance$gene

length(unique(sample_performance$gene))
sample_performance_merge$Model = str_replace(sample_performance_merge$mode, 'random', 'Random')

concise_sample_performance = sample_performance_merge[,c('rsq', 'mode', 'selected_features')]
print(concise_sample_performance)


table(sample_performance_merge$mode)

four_modes_plot <-ggplot(sample_performance_merge, aes(x = performance , color = mode )) + geom_density()
ggplot(subset(sample_performance_merge, mode == 'AlltfShuffle' & rsq > 0), aes(x = rsq , color = mode )) + geom_density()


TF_random_plot <- ggplot(subset(sample_performance_merge, mode != 'All' & mode != 'SNP'),
                         aes(x = performance, color = mode )) + geom_density()
ggsave(f_p('%s/TF_random_performance.tiff', results_dir), plot = TF_random_plot)

sample_performance_merge


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
all_feature_table = f_summary_selected_features(loc_batch_name, c('chr22'), mode_list[3], output_feature_file = 'feature.table', threshold )

SNP_feature_table = f_summary_selected_features(loc_batch_name, c('chr22'), mode_list[5], output_feature_file = 'feature.table.tf', threshold )



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

################Best genes#####################
best_index=which.max(sample_performance_merge$performance)
best_gene = sample_performance_merge[best_index,'gene']
sample_performance_merge %>% filter(gene == best_gene)




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
