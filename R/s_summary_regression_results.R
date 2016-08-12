
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
    make_option(c("--collection"),      type="character", default='addPenalty', help="The parameter based name", metavar="character"),
    make_option(c("--chr_batch"),      type="character", default='1chr', help="1chr for chr22, 3chr for (2, 10, 22)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$batch_name)){
    #batch_name = '54samples_evalue'
    #batch_name = '462samples_sailfish_quantile'
    #batch_name = '462samples_snyder_norm'
    #batch_name = '462samples_quantile'
    #batch_name = '462samples_quantile_rmNA'
    #batch_name = '462samples_log_quantile'
    batch_name = '445samples_sailfish'
    ##chr_str = 'chr22'
    chr_str = 'chr22'
}else{
    batch_name = opt$batch_name
    chr_str = opt$chr
    colletion_name = opt$collection
}

modes_list = list()

modes_list$non_pop = c(All='rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.all',
                 SNP='rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.SNP',
                 TF='rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.TF',
                 random = 'rm.histone_rm.miRNA_model.glmnet_rm.predict.TF_add.YRI_population.all_TF.exp.type.fakeTF_add.gm12878_new.batch.462samples.snyder.original_batch.mode.random')


modes_list$pop = c(
    AlltfShuffle ='rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.AlltfShuffle',
    noInteract = 'rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.noInteract',
    All             = 'rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.all',
             SNP = 'rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.SNP',
             SNPinTF = 'rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.SNPinTF',
    fakeInteract = 'rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.fakeInteract',
             #random = 'rm.histone_rm.miRNA_model.glmnet_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.random',
             TF = 'rm.histone_rm.miRNA_model.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.TF')

modes_list$cvglmnet = c(
    AlltfShuffle ='rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.AlltfShuffle',
    #noInteract = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.noInteract',
    All             = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.all',
    #AllnoInteract = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.AllnoInteract',
    SNP = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.SNP',
    SNPinTF = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.SNPinTF',
    fakeInteract = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.fakeInteract',
             #random = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.random',
             TF = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.TF',
             TFShuffle = 'rm.histone_rm.miRNA_model.cv.glmnet_rm.penalty_population.all_add.gm12878_new.batch.462samples.snyder.original_batch.mode.TFShuffle'
)

modes_list$rmCor=c(
All            = 'rm.histone_model.cv.glmnet_rm.penalty_population.all_new.batch.462samples.snyder.original_batch.mode.All_other.info.rmCor',
SNP         = 'rm.histone_model.cv.glmnet_rm.penalty_population.all_new.batch.462samples.snyder.original_batch.mode.SNP_other.info.rmCor',
SNPinTF  = 'rm.histone_model.cv.glmnet_rm.penalty_population.all_new.batch.462samples.snyder.original_batch.mode.SNPinTF_other.info.rmCor'
)

modes_list$keepCor=c(
All = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.All_other.info.normCor',
SNPinTF = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.SNPinTF_other.info.normCor',
SNP = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.SNP_other.info.normCor',
InterOnlySNPinTF='rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.InterOnlySNPinTF_other.info.normCor'
)

modes_list$rmPenalty=c(
All = 'rm.histone_model.cv.glmnet_rm.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.All_other.info.normCor',
SNPinTF = 'rm.histone_model.cv.glmnet_rm.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.SNPinTF_other.info.normCor',
SNP = 'rm.histone_model.cv.glmnet_rm.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.SNP_other.info.normCor',
TF = 'rm.histone_model.cv.glmnet_rm.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.TF_other.info.normCor',
InterOnlySNPinTF='rm.histone_model.cv.glmnet_rm.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.InterOnlySNPinTF_other.info.normCor'
)


modes_list$addPenalty=c(
All = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.All_other.info.normCor',
#AllfilterMinor = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.AllfilterMinor_other.info.normCor',
#AlltopTF = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.AlltopTF_other.info.normCor',
#All = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.All_other.info.normCor',
SNPinTF = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.SNPinTF_other.info.normCor',
SNP = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.SNP_other.info.normCor',
TF = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.TF_other.info.normCor'
,
#AllrmHic='rm.histone_model.cv.glmnet_rm.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.All_other.info.normCor',
#AllnoInteract = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.AllnoInteract_other.info.normCor',
##InterOnlySNPinTF='rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.InterOnlySNPinTF_other.info.normCor',
#noInteract ='rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.noInteract_other.info.normCor',
AlltfShuffle = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.AlltfShuffle_other.info.normCor'
)

modes_list$snpOnly=modes_list$addPenalty[c('All', 'TF', 'SNP')]

#loc_batch_name = '462samples_quantile_rmNA'
#collection_name = 'rmCor'



modes_list$addPenaltyRmCor = f_create_new_mode_list(modes_list$addPenalty, 'normCor', 'rmCor')
modes_list$snyder.norm = f_create_new_mode_list(modes_list$addPenalty, 'snyder.original', 'snyder.norm')
modes_list$maxit = f_create_new_mode_list(modes_list$addPenalty, 'normCor', 'maxit1M')
modes_list$lm = f_create_new_mode_list(modes_list$addPenalty, 'normCor', 'lm')

loc_batch_name = '445samples_sailfish'
#loc_batch_name = '445samples_snpOnly'
collection_name='snyder.norm'
#mode_list = modes_list[[collection_name]][c(1,3)]
mode_list = modes_list[[collection_name]][c('All', 'TF', 'SNP')]
options(width=60)
#' #Modes used in the analysis
print(str_replace_all(mode_list, 'rm.histone_model.|other|new.batch', ''), width = 50)



#loc_batch_name = '462samples_genebody'
#mode_list = genebody_list
#chr_num_list = c(22)
chr_num_list = c(22, 2, 10)

f_built_file_name_from_paras('s_summary_regression_results', f_p('%s_%s_%schrs', loc_batch_name, collection_name, length(chr_num_list)))



#' #Collect the performance data.
sample_performance_merge = f_collect_performance_in_multiple_chrs(loc_batch_name, mode_list, chr_num_list)
sample_performance_merge$performance = sample_performance_merge$rsq
results_dir = f_p('./data/%s/rnaseq/results/', loc_batch_name)



#' #Thresholds
table(sample_performance_merge$mode)
#threshold = quantile( subset(sample_performance_merge, mode == 'random')$performance , probs = 0.99 )
threshold = 0.05



#' #Stats For the robust genes
sample_performance_merge %>% filter( performance > 0.8 * train_performance, performance > 0.05 ) %>% group_by(mode) %>%
    dplyr::summarise(good_predictions = sum(performance > threshold),
                               mean_performance = mean(performance),
                               mean_SD = mean(SD),
                               feature_num = mean(selected_features),
                     max_feature = max(selected_features))

head(sample_performance_merge)

table(sample_performance_merge$mode)


#' #Stats For the intersect genes
gene_count = table(sample_performance_merge$gene) == length(mode_list)
intersect_genes = names(gene_count)[gene_count]

intersect_stats <- sample_performance_merge %>% filter(performance > threshold, gene %in% intersect_genes) %>% group_by(mode) %>%
    dplyr::summarise(Total = length(intersect_genes),
                               good_predictions = sum(performance > threshold),
                               mean_performance = mean(performance),
                               mean_train = mean(train_performance),
                               mean_SD = mean(SD),
                               mean_input_feature = mean(num_feature), 
                               mean_selected_features = mean(selected_features),
                               max_feature = max(selected_features))
print(intersect_stats)
#write.table(format(as.data.frame(intersect_stats), digits = 3), f_p('%s/basic_stats.txt', results_dir), quote = FALSE, sep = '\t', row.names = FALSE)


#' Ensemble all the models
##sample_performance_merge %>% filter(performance > 0.6) %>% arrange(gene)
target_modes = c('All', 'TF', 'SNP', 'SNPinTF')
ensemble_mode=sample_performance_merge %>% filter(gene %in% intersect_genes, mode %in% target_modes,
                                                  performance > threshold) %>% group_by(gene) %>%
    dplyr::summarise(performance = max(rsq), index = mode[which.max(rsq)])
cat('Distribution of selected models:\n')
table(ensemble_mode$index)
snp_mode = sample_performance_merge %>% filter(mode == 'SNP', gene %in% ensemble_mode$gene)
mean_improvement=(mean(ensemble_mode$performance) - mean(snp_mode$performance))/mean(snp_mode$performance)
cat('Mean improvment over the SNP model', mean_improvement, '\n')
if ('AllrmHic' %in% sample_performance_merge$mode){
    rmHic = subset(sample_performance_merge, mode == 'AllrmHic' & performance > 0.05)
    rownames(rmHic) = rmHic$gene
    dim(rmHic)
    all = subset(sample_performance_merge, mode == 'All')
    rownames(all) = all$gene
    mean_improvement=(mean(all[rownames(rmHic),'performance'], na.rm = TRUE) - mean(rmHic$performance))/mean(rmHic$performance)
    cat('Mean improvment over the SNP model', mean_improvement, '\n')
}

#sample_performance_merge %>% filter(gene == as.data.frame(ensemble_mode)[2,'gene'])




print(intersect_stats)

#' ##Plots to compare the performance
f_venn_plot_overlap(sample_performance_merge, intersect_genes, threshold)

table(sample_performance_merge$mode)

good_performance <- sample_performance_merge %>% filter(gene %in% intersect_genes, rsq > 0.05)
dim(sample_performance_merge)

head(good_performance)

good_performance2 = good_performance %>% filter(mode == 'SNP' | mode == 'All', performance > 0.05)
ggplot(good_performance2, aes(x = performance , color = mode )) + geom_density()

head(good_performance2)
group_A  = 'AlltfShuffle'
group_B  = 'All'
thres = 0.005





f_compare_improvment_for_two_groups('SNP', "All", good_performance)
f_compare_improvment_for_two_groups('SNP', "TF", good_performance)
f_compare_improvment_for_two_groups('SNP', "SNPinTF", good_performance)
f_compare_improvment_for_two_groups('SNPinTF', "TF", good_performance)


f_plot_performance_and_stats_test(good_performance, 'SNP', 'All' )
f_plot_performance_and_stats_test(good_performance, 'AllrmHic', 'All')
f_plot_performance_and_stats_test(good_performance, 'SNPinTF', 'TF' )
f_plot_performance_and_stats_test(good_performance, 'SNP', 'TF' )
f_plot_performance_and_stats_test(good_performance, 'All', 'TF' )
#f_plot_performance_and_stats_test(good_performance, 'SNP', 'AllnoInteract' )
f_plot_performance_and_stats_test(good_performance, 'SNPinTF', 'SNP' )
f_plot_performance_and_stats_test(good_performance, 'SNPinTF', 'TF' )
f_plot_performance_and_stats_test(good_performance, 'AllnoInteract', 'All' )
f_plot_performance_and_stats_test(good_performance, 'SNP', 'AlltfShuffle')
f_plot_performance_and_stats_test(good_performance, 'AlltfShuffle', 'All' )
f_plot_performance_and_stats_test(good_performance, 'TF', 'fakeInteract' )
f_plot_performance_and_stats_test(good_performance, 'TF', 'TFShuffle')
f_plot_performance_and_stats_test(good_performance, 'TF', 'noInteract')
#f_plot_performance_and_stats_test(good_performance, 'All', 'TF')
#f_plot_performance_and_stats_test(good_performance, 'TF', 'noInteract')
#f_plot_performance_and_stats_test(good_performance, 'TF', 'fakeInteract')
#f_plot_performance_and_stats_test(good_performance, 'noInteract', 'fakeInteract') #Not significant


#' #Test performance drop against model variance
try(f_test_performance_drop_with_ml_variance(sample_performance_merge, 'SNP', 'SNPinTF'))

sample_performance_merge$Model = str_replace(sample_performance_merge$mode, 'random', 'Random')
ggplot(subset(sample_performance_merge, performance > 0), aes(x = performance , color = mode )) + geom_density()


#' #GO analysis, wait until the whole genome analysis.
source('s_GO_analysis_predictions.R')
#thres is working.
f_performance_go_analysis(sample_performance_merge, 'TF', thres = 0.2) #check this after the whole genome analysis.


