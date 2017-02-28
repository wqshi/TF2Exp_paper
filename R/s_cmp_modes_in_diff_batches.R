#+ fig.height=8

source('s_summary_fun.R')
library(gridExtra)
library(tidyr)
source('~/R/s_ggplot2_theme.R', chdir = TRUE) 


chr_mode = 'chr22'
normal_batch = '358samples_regionkeepLow'
rsync_flag = TRUE


if (chr_mode == 'chr22'){
    chr_list = c('chr22')
}else if (chr_mode == '3chrs' ){
    chr_list = c('chr22', 'chr10', 'chr15')
}




f_compare_two_modes_in_diff_batches <- function(batch_A, batch_B, mode_index1, mode_index2, modes_list1, modes_list2, chr_list, perf_thres = 0.05, rsync_flag = F, performance_col = 'performance', add_batch_name = TRUE, debug = FALSE){

    if (f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
        browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
    }
    
    performance_merge1 = data.frame()
    performance_merge2 = data.frame()
    for (chr_str in chr_list){
        dir.create(f_p('./data/%s/rnaseq/%s/', batch_A, chr_str))
        return1 =f_summary_regression_results(batch_A, chr_str = chr_str, mode_name = modes_list1[[mode_index1]], return_features = T, rsync_flag = rsync_flag)
        performance_merge1 = rbind(performance_merge1, return1$performance)

        dir.create(f_p('./data/%s/rnaseq/%s/', batch_B, chr_str))
        return2 =f_summary_regression_results(batch_B, chr_str = chr_str, mode_name = modes_list2[[mode_index2]], return_features = T, rsync_flag = rsync_flag)
        performance_merge2 = rbind(performance_merge2, return2$performance)
    }

    if (add_batch_name){
        loc_batch_A = paste0('A-', batch_A, ':',  (mode_index1))
        loc_batch_B = paste0('B-', batch_B, ':',  (mode_index2))
    }else{
        loc_batch_A = (mode_index1 )
        loc_batch_B = (mode_index2)
    }

    performance_merge1$mode = loc_batch_A
    performance_merge2$mode = loc_batch_B
    
    shared_cols=intersect(colnames(performance_merge1), colnames(performance_merge2))

    performance_merge = rbind(performance_merge1[, shared_cols], performance_merge2[, shared_cols])
    performance_merge$performance = performance_merge[,performance_col]
    
    dim(performance_merge2)
    dim(performance_merge1)

    setdiff(rownames(performance_merge1), rownames(performance_merge2) )
    setdiff(rownames(performance_merge2), rownames(performance_merge1) )

    table(performance_merge$mode)
    dot_plot=f_plot_performance_and_stats_test(subset(performance_merge, performance > perf_thres), loc_batch_A, loc_batch_B)
    compare_results = f_compare_improvment_for_two_groups(loc_batch_B, loc_batch_A, subset(performance_merge, performance > perf_thres), thres = 0.01, return_flag = F, debug = F)

    #final_plot <-arrangeGrob( arrangeGrob(dot_plot, compare_results$diff, nrow =1),  compare_results$overall_density, compare_results$density, nrow = 3)
    final_plot <-arrangeGrob( arrangeGrob(dot_plot, compare_results$diff, nrow =1), compare_results$density, nrow = 2)
    plot(final_plot)

    return (performance_merge)
}





#print(str_replace_all(keepZero_list, 'rm.histone_model.|other|new.batch.', ''), width = 50)


#normal_batch = '445samples_region'
#normal_batch = '445samples_diff'
snp_batch = f_p('%s_snpOnly', normal_batch)








##################Re################
chr_list = c('chr22', 'chr10', 'chr15')
#chr_list = c('chr22')



collection_name = 'peer358cor'
peer_list = modes_list[[collection_name]][c('TF', 'SNP','All', 'SNPinTF', 'random' ,'TFaddInteract','TFsnpMatch', 'TFaddPenalty', 'TFfilterMinor')]
rmdup_list = modes_list[['peer358corRmdup']][c('TF', 'SNP','All', 'SNPinTF', 'random' ,'TFaddInteract','TFsnpMatch', 'TFaddPenalty', 'TFfilterMinor')]
gtex_list = modes_list[['gtex']][c('TF', 'SNP','All', 'SNPinTF', 'random' ,'TFaddInteract','TFsnpMatch', 'TFaddPenalty', 'TFfilterMinor')]
elastic_list = modes_list[['elastic']][c('TF', 'SNP','All', 'SNPinTF', 'random' ,'TFaddInteract','TFsnpMatch', 'TFaddPenalty', 'TFfilterMinor')]

#Compare the performance between lasso and elastic
tf_performance = f_compare_two_modes_in_diff_batches(normal_batch, normal_batch, 'TF', 'TF', rmdup_list, elastic_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)


tf_performance = f_compare_two_modes_in_diff_batches(normal_batch, snp_batch, 'TFsnpMatch', 'TF', gtex_list, gtex_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)




tf_performance = f_compare_two_modes_in_diff_batches(normal_batch, normal_batch, 'TF', 'TF', peer_list, rmdup_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)

tf_performance = f_compare_two_modes_in_diff_batches(normal_batch, normal_batch, 'TF', 'TFfilterMinor', rmdup_list, rmdup_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)

tf_performance = f_compare_two_modes_in_diff_batches(normal_batch, normal_batch, 'TF', 'TF', rmdup_list, gtex_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)

tf_performance = f_compare_two_modes_in_diff_batches(normal_batch, normal_batch, 'TF', 'All', gtex_list, gtex_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)

tf_performance = f_compare_two_modes_in_diff_batches(normal_batch, normal_batch, 'TF', 'All', rmdup_list, rmdup_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)

tf_match_performance = f_compare_two_modes_in_diff_batches(normal_batch, snp_batch, 'TFsnpMatch', 'TF', rmdup_list, rmdup_list, chr_list, rsync_flag = rsync_flag, add_batch_name = T)


t_test_feature_counts <- function(tf_match_performance){
    feature_counts <- tf_performance %>% select(gene, mode, num_feature) %>% spread(mode, num_feature) %>% head
    f_ASSERT(all(feature_counts[,2] == feature_counts[,3]))
}

t_test_feature_counts(tf_match_performance)


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

