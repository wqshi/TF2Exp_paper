#+ fig.height=8

source('s_summary_fun.R')
library(gridExtra)
library(tidyr)
mode_list = c(TF = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.TF_other.info.tradR2',
              SNP = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.SNP_other.info.tradR2',
              SNPinTF = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.SNPinTF_other.info.tradR2',
              All = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.All_other.info.tradR2',
              random = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.random_other.info.tradR2',
              noInteract = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.noInteract_other.info.tradR2',
              TFsnpMatch = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.TFsnpMatch_other.info.tradR2')


keepZero_list = f_create_new_mode_list(mode_list, 'add.penalty', 'rm.penalty_rm.YRI')

source('~/R/s_ggplot2_theme.R', chdir = TRUE) 



chr_mode = '3chrs'
normal_batch = '358samples_regionkeepLow'
pop = 'Pop' #or false
rsync_flag = TRUE


if (chr_mode == 'chr22'){
    chr_list = c('chr22')
}else if (chr_mode == '3chrs' ){
    chr_list = c('chr22', 'chr10', 'chr15')
}

if (pop == 'Pop'){
    keepZero_list = f_create_new_mode_list(keepZero_list, 'tradR2', 'tradR2keepZeroPop')
    keepZero_list = f_create_new_mode_list(keepZero_list, '445', '358')
}else{
    keepZero_list = f_create_new_mode_list(keepZero_list, 'tradR2', 'tradR2keepZero')
}

print(str_replace_all(keepZero_list, 'rm.histone_model.|other|new.batch.', ''), width = 50)

f_built_file_name_from_paras('s_cmp_modes_in_diff_batches', f_p('%s_%s_%schrs', normal_batch, pop, length(chr_list)))
#normal_batch = '445samples_region'
#normal_batch = '445samples_diff'
snp_batch = f_p('%s_snpOnly', normal_batch)





##################Re

tf_performance = f_compare_modes_in_diff_batches(normal_batch, snp_batch, 'TFsnpMatch', 'TF', keepZero_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)
library(tidyr)
feature_counts <- tf_performance %>% select(gene, mode, num_feature) %>%  spread(key = mode, value = num_feature )
f_ASSERT(all(feature_counts[2] <= feature_counts[3]), 'Feature nums are different in 445samples_diff and 445samples_diff_snp_only')
head(feature_counts)

tf_performance = f_compare_modes_in_diff_batches(normal_batch, normal_batch, 'TFsnpMatch', 'TF', keepZero_list, chr_list, rsync_flag = rsync_flag)
tf_performance = f_compare_modes_in_diff_batches(normal_batch, normal_batch, 'random', 'TF', keepZero_list, chr_list, rsync_flag = TRUE, perf_thres = -1)

#tf_performance = f_compare_modes_in_diff_batches(normal_batch, snp_batch, 'TF', 'TF', keepZero_list, chr_str, rsync_flag = T)
tf_performance = f_compare_modes_in_diff_batches(normal_batch, normal_batch, 'All', 'SNP', keepZero_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)
tf_performance = f_compare_modes_in_diff_batches(normal_batch, normal_batch, 'TF', 'SNP', keepZero_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)
tf_performance = f_compare_modes_in_diff_batches(normal_batch, normal_batch, 'SNPinTF', 'SNP', keepZero_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)
tf_performance = f_compare_modes_in_diff_batches(normal_batch, normal_batch, 'TFsnpMatch', 'SNPinTF', keepZero_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)
tf_performance = f_compare_modes_in_diff_batches(normal_batch, normal_batch, 'TF', 'SNPinTF', keepZero_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)
#tf_performance = f_compare_modes_in_diff_batches(normal_batch, normal_batch, 'TF', 'SNPinTF', keepZero_list, chr_list, rsync_flag = T)
#tf_performance = f_compare_modes_in_diff_batches(normal_batch, normal_batch, 'TF', 'SNP', keepZero_list, chr_list, rsync_flag = T, add_batch_name =F, debug = F)

table(tf_performance$mode)
colnames(tf_performance)









if (FALSE){


    feature_counts %>% filter(gene == 'ENSG00000100124.8')
    feature_counts %>% filter(gene == 'ENSG00000128274.11')
    sum(feature_counts[,3] >= feature_counts[,2])/nrow(feature_counts)
    which(feature_counts[,3] - feature_counts[,2] < 0)
    
    
    f_compare_modes_in_diff_batches('445samples_diff', '445samples_diff_snpOnly', 'TFsnpMatch', 'TF' , keepZero_list, chr_str, rsync_flag = T)
    
    f_compare_modes_in_diff_batches('445samples_diff', '445samples_diff', 'SNP', 'All', keepZero_list, chr_str, rsync_flag = T)
    
    f_compare_modes_in_diff_batches('445samples_maxVar', '445samples_maxVar', 'TF', 'TFsnpMatch', mode_list, chr_str, rsync_flag = T)
    
    f_compare_modes_in_diff_batches('445samples_sailfish', '445samples_diff', 'TFsnpMatch', 'TFsnpMatch', mode_list, chr_str)
    
    f_compare_modes_in_diff_batches('445samples_diff', '445samples_sailfish', 'TF', 'TF', mode_list, chr_str)
}





if (FALSE){

    f_compare_modes_in_diff_batches(normal_batch, snp_batch, 'TF', 'TF', keepZero_list, chr_str, rsync_flag = T)
    f_compare_modes_in_diff_batches(normal_batch, normal_batch, 'All', 'SNP', keepZero_list, chr_str, rsync_flag = T)
    a = f_compare_modes_in_diff_batches(normal_batch, snp_batch, 'TFsnpMatch', 'TF', keepZero_list, chr_str, rsync_flag = T)
    f_compare_modes_in_diff_batches(normal_batch, normal_batch, 'TF', 'SNPinTF', keepZero_list, chr_str, perf_thres = 0.05, rsync_flag = T)

    snp_performance = f_compare_modes_in_diff_batches(normal_batch, snp_batch, 'SNPinTF', 'SNPinTF', keepZero_list, chr_str, perf_thres = 0.05, rsync_flag = T)

    head(snp_performance)

    snp_performance[which.max(snp_performance$performance),]
    snp_performance[which.min(snp_performance$performance),]


    snp_performance = f_compare_modes_in_diff_batches('445samples_diff', '445samples_diff_snpOnly', 'SNP', 'SNP', keepZero_list, chr_str, rsync_flag = T)
    perf_diff <- snp_performance %>% select(gene, mode, performance) %>% spread(mode, performance)
    perf_diff$diff = perf_diff[,2] - perf_diff[,3]
    perf_diff = f_sort_by_col(perf_diff, 'diff', decr_flag = F)

    head(perf_diff)
    tail(perf_diff)

    table(snp_performance$mode)

    TF_performance = f_compare_modes_in_diff_batches(normal_batch, snp_batch, 'TFsnpMatch', 'TF', keepZero_list, chr_str, perf_thres = 0.05 , rsync_flag = T)

    perf_diff <- tf_performance %>% select(gene, mode, performance) %>% spread(mode, performance)

    perf_diff$diff = perf_diff[,2] - perf_diff[,3]

    perf_diff = f_sort_by_col(perf_diff, 'diff', decr_flag = T)

    sum(perf_diff$diff < 0, na.rm = T)
    sum(perf_diff$diff > 0, na.rm = T)

    head(perf_diff)

    TF_performance %>% filter(gene == as.character(perf_diff$gene[1]))

    max_gene = rownames(compare_results[which.max(compare_results$perf_diff),])

    return2$features[grep(max_gene, return2$features$name),]
    return1$features[grep(max_gene, return1$features$name),]

    ##Get the min and max performance genes.
    return1 =f_summary_regression_results('445samples_region', chr_str = chr_str, mode_name = keepZero_list['SNPinTF'], return_features = T)
    max_gene = return1$performance[which.max(return1$performance$performance), ]
    min_gene = return1$performance[which.min(return1$performance$performance), ]
}
