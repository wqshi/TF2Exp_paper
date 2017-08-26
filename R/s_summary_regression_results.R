##This script combine the results from s_regression_for_one_gene.R
#Sync the data from clustdell.
library("optparse")
option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="like 445samples_sailfish, 445samples_rareVar", metavar="character"),
    make_option(c("--collection"),      type="character", default='addPenalty', help="e.g. addPenalty, snyder.norm", metavar="character"),    
    make_option(c("--chr_batch"),      type="character", default='1chr', help="1chr for chr22, 3chrs for (2, 10, 15)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
str(opt)


library(reshape2)
library(futile.logger)
library(dplyr)
library(matrixStats)
setwd('~/expression_var/R/')
source('~/R/s_function.R', chdir = TRUE)
suppressMessages(source('s_gene_regression_fun.R'))
source('s_summary_fun.R')
library(stringr)
library(gridExtra)
source('~/R/s_ggplot2_theme.R', chdir = TRUE)
source('./s_project_data.R')
source('s_project_funcs.R')
options(dplyr.width = Inf)



if (is.null(opt$batch_name)){
    #batch_name = '54samples_evalue'
    #batch_name = '462samples_sailfish_quantile'
    #batch_name = '462samples_snyder_norm'
    #batch_name = '462samples_quantile'
    #batch_name = '462samples_quantile_rmNA'
    #batch_name = '462samples_log_quantile'
    #loc_batch_name = '445samples_rareVar'
    loc_batch_name = '358samples_regionkeepLow'
    #loc_batch_name = '358samples_regionkeepLow'
    #loc_batch_name = '445samples_maf'
    ##chr_str = 'chr22'
    chr_num_list = c(1:22)
    #chr_num_list = c(22)
    #collection_name = 'non_pop'
    #collection_name = 'peer358cor'
    #collection_name = 'rmZero'
    collection_name = 'peer358corRmdup'    
}else{
    loc_batch_name = opt$batch_name

    if (opt$chr_batch == '1chr'){
        chr_num_list = c(5)
    }else{
        chr_num_list = 1:22
    }
    collection_name = opt$collection
    
}

R2_type = 'cor' # 'tradR2'
print('Chr_num_list:')
print(chr_num_list)
print(collection_name)

batch_size = str_replace(loc_batch_name, 'samples.*', '')
#mode_list defined in the s_summary_fun22	16051249	rs62224609	T	C
mode_list = modes_list[[collection_name]][c('TF', 'SNP','All', 'SNPinTF', 'random' ,'TFaddInteract','TFsnpMatch', 'TFaddPenalty')]
#mode_list = modes_list[[collection_name]][c('TF', 'SNP')]

mode_list
#' #Modes used in the analysis
print(str_replace_all(mode_list, 'rm.histone_model.|other|new.batch.', ''), width = 50)

f_built_file_name_from_paras('s_summary_regression_results', f_p('%s_%s_%s_%schrs', loc_batch_name, collection_name, R2_type ,length(chr_num_list)))
results_dir = f_p('./data/%s/rnaseq/results/%s/', loc_batch_name, collection_name)
figure_dir = f_p('./data/%s/rnaseq/results/%s/figures', loc_batch_name, collection_name)
dir.create(results_dir)


perf_file = f_p('%s/performance_results.txt', results_dir)

source('s_test_cat.R')
f_write_paper_results('=====Performance file======', data = date(), file = perf_file)

cat(collection_name, loc_batch_name, '\n', file = perf_file, append = T)

#####################Collect the performance data.###############
read_data = F
if (read_data == TRUE){

    sample_performance_merge = f_collect_performance_in_multiple_chrs(loc_batch_name, mode_list, chr_num_list)

    snpOnly_performance_merge = f_collect_performance_in_multiple_chrs('358samples_regionkeepLow_snpOnly', mode_list['TF'], chr_num_list)
    snpOnly_performance_merge$mode = paste0( snpOnly_performance_merge$mode, '.snp')
    sample_performance_merge = rbind(sample_performance_merge, snpOnly_performance_merge)
    
    sample_performance_merge$performance = sample_performance_merge$rsq
    write.table(sample_performance_merge, file = f_p('%s/sample_performance_merge', results_dir), quote = F, sep = '\t', row.names = F)

    #Write the max performance
    #max_perf_df = sample_performance_merge %>% group_by(gene) %>% summarise(max_perf = max(performance)) %>% filter(max_perf > 0)
    #flog.info( "%s out of %s genes have performance abover 0",  sum(max_perf_df$max_perf > 0), nrow(max_perf_df) )
    #write.table(max_perf_df, file = './data/r_results/max_perf.txt', quote = F, row.names = F)
    
}else{
    sample_performance_merge = read.table(file = f_p('%s/sample_performance_merge', results_dir), header = T)
}

######################Collect data done###############################

head(sample_performance_merge)





######################1. Overall stats of the performance in different modes#########################################
#' #Thresholds
table(sample_performance_merge$mode, sample_performance_merge$chr)
threshold = 0.05 #Given by other paper


#' #Stats For the intersect genes
#sample_performance_merge <- sample_performance_merge %>% filter(!mode %in% c('random', 'TFaddInteract'))

gene_count = table(sample_performance_merge$gene) == length(unique(sample_performance_merge$mode))
intersect_genes = names(gene_count)[gene_count]
rest_genes = names(gene_count)[!gene_count]

if (FALSE){
    ##For the manual check
    sample_performance_merge %>% filter(gene %in% rest_genes)
    sample_performance_merge %>% filter(gene == 'ENSG00000073169.9')
    setdiff( (sample_performance_merge %>% filter( mode == 'TFsnpMatch'))$gene, (sample_performance_merge %>% filter( mode == 'SNP'))$gene )
    setdiff( (sample_performance_merge %>% filter( mode == 'SNP'))$gene, (sample_performance_merge %>% filter( mode == 'TF'))$gene )
}
      
intersect_stats <- sample_performance_merge %>% filter(performance > threshold, gene %in% intersect_genes) %>% group_by(mode) %>%
    dplyr::summarise(Total = length(intersect_genes),
                     good_predictions = sum(performance > threshold),
                     mean_performance = mean(performance),
                     mean_train = mean(train_performance),
                     mean_SD = mean(SD),
                     mean_input_feature = mean(num_feature), 
                     mean_selected_features = mean(selected_features),
                     max_feature = max(selected_features),
                     CHB_count = sum(!is.na(CHB)), CHB_mean = mean(CHB, na.rm = T),
                     JPT_count = sum(!is.na(JPT)), JPT_mean = mean(JPT, na.rm = T)
                     )

#Output: mean validation in microarray.TF and SNPinTF. In theory should be TF.snp.
intersect_stats

basic_stats <- sample_performance_merge %>% group_by(mode) %>%
    dplyr::summarise(Total = length(mode),
                     good_predictions = sum(performance > threshold),
                     mean_performance = mean(performance),
                     mean_train = mean(train_performance),
                     mean_SD = mean(SD),
                     mean_input_feature = mean(num_feature), 
                     mean_selected_features = mean(selected_features),
                     max_feature = max(selected_features),
                     CHB_count = sum(!is.na(CHB)), CHB_mean = mean(CHB, na.rm = T),
                     JPT_count = sum(!is.na(JPT)), JPT_mean = mean(JPT, na.rm = T)
                     )

print(intersect_stats)
tf_stats <- sample_performance_merge %>% filter(mode == 'TF') %>%
    dplyr::summarise(Total = length(mode),
                     good_predictions = sum(performance > threshold),
                     good_ratio = sum(performance > threshold)/ length(mode),
                     mean_performance = mean(performance),
                     mean_train = mean(train_performance),                     
                     mean_SD = mean(SD),
                     mean_hic = mean(hic_num),
                     mean_input_feature = mean(num_feature), 
                     mean_selected_features = mean(selected_features, na.rm = T),
                     max_feature = max(selected_features, na.rm = T),
                     CHB_count = sum(!is.na(CHB)), CHB_mean = mean(CHB, na.rm = T),
                     JPT_count = sum(!is.na(JPT)), JPT_mean = mean(JPT, na.rm = T)
                     )
tf_stats

tf_good_stats <- sample_performance_merge %>% filter(mode == 'TF', performance > threshold) %>%
    dplyr::summarise(Total = length(mode),
                     mean_performance = mean(performance),
                     mean_train = mean(train_performance),                     
                     mean_SD = mean(SD),
                     mean_input_feature = mean(num_feature),
                     total_selected_features = sum(selected_features),
                     mean_selected_features = mean(selected_features, na.rm = T),
                     max_feature = max(selected_features, na.rm = T),
                     CHB_count = sum(!is.na(CHB)), CHB_mean = mean(CHB, na.rm = T),
                     JPT_count = sum(!is.na(JPT)), JPT_mean = mean(JPT, na.rm = T)
                     )

tf_good_stats

random_stats <- sample_performance_merge %>% filter(mode == 'random') %>%
    dplyr::summarise(Total = length(mode),
                     good_predictions = sum(performance > threshold),
                     good_ratio = sum(performance > threshold)/ length(mode),
                     mean_performance = mean(performance),
                     max_performance = max(performance)
                     )


#sample_performance_merge %>% filter(performance > threshold, mode == 'TF' ,gene %in% intersect_genes) %>% select(gene, performance, CEU, CHB, JPT) 

table(sample_performance_merge$mode)

#print(intersect_stats)
f_write_paper_results('Raw numbers from TF models', tf_stats, perf_file, scientific = T)
f_write_paper_results('Paper numbers',
                      f_p('Total %s, good %s, average hic %s, average TF binding events %s',
                          tf_stats$Total, tf_stats$good_predictions, tf_stats$mean_hic, tf_stats$mean_input_feature), perf_file)

write.table(format(as.data.frame(basic_stats), digits = 3), f_p('%s/basic_stats.txt', results_dir), quote = FALSE, sep = '\t', row.names = FALSE)
write.table(format(as.data.frame(tf_good_stats), digits = 3), f_p('%s/predictable_tf_stats.txt', results_dir), quote = FALSE, sep = '\t', row.names = FALSE)
write.table(format(as.data.frame(random_stats), digits = 3), f_p('%s/random_stats.txt', results_dir), quote = FALSE, sep = '\t', row.names = FALSE)
write.table(format(as.data.frame(tf_good_stats), digits = 3), f_p('%s/tf_good_stats.txt', results_dir), quote = FALSE, sep = '\t', row.names = FALSE)

##########################Overall stats done#########################################







########################Section 2: Ensemble all the models###########################
##sample_performance_merge %>% filter(performance > 0.6) %>% arrange(gene)
target_modes = c('All', 'TF', 'SNP', 'SNPinTF', 'random')
ensemble_mode=sample_performance_merge %>% filter(gene %in% intersect_genes, mode %in% target_modes,
                                                  performance > threshold) %>% group_by(gene) %>%
    dplyr::summarise(performance = max(rsq), index = mode[which.max(rsq)])
cat('Distribution of selected models:\n')
table(ensemble_mode$index)
snp_mode = sample_performance_merge %>% filter(mode == 'SNP', gene %in% ensemble_mode$gene)
mean_improvement=(mean(ensemble_mode$performance) - mean(snp_mode$performance))/mean(snp_mode$performance)
cat('Mean improvment over the SNP model', mean_improvement, '\n')
##########################Ensemble done###############################################





##############Section 3: Compare adding HiC and without HiC###############
if ('TFaddPenalty' %in% sample_performance_merge$mode){
    rmHic = subset(sample_performance_merge, mode == 'TFaddPenalty' & performance > 0.05)
    rownames(rmHic) = rmHic$gene
    dim(rmHic)
    all = subset(sample_performance_merge, mode == 'TF')
    rownames(all) = all$gene
    mean(all[rownames(rmHic),'performance'] - rmHic$performance)

    mean(rmHic$performance)
    cat('Mean improvment over the SNP model', mean_improvement, '\n')
}
##################Compare adding Hic done########################################







####################Section 4. Pairwise comparison of two mdoels#######################
#f_venn_plot_overlap(sample_performance_merge, intersect_genes, threshold)

good_performance <- sample_performance_merge %>% filter(gene %in% intersect_genes, rsq > 0.05)

library(scales)
#Random vs TF plots
performance_sort <- sample_performance_merge %>% arrange(mode, performance) %>%
    group_by(mode) %>% dplyr::mutate( rank = row_number(mode)/length(mode) ) %>%
    filter(mode %in% c('TF', 'random', 'TFaddPenalty', 'TFaddInteract'))
random_tf_plot <-
    ggplot(performance_sort, aes(rank, performance, color = mode)) + geom_point(size = 0.1) +
    theme_Publication( base_size = 15) + 
    theme(legend.position = c(0.2,0.8), legend.direction = 'vertical') + geom_hline(yintercept=0.05, linetype = '1F') +
    xlab('Cumulative percentages of investigated genes') + ylab('Model Performance') + 
    #scale_y_continuous(breaks=c(0,0.05,0.2,0.4,0.6,0.8)) +
    scale_x_continuous(labels=percent, limits = c(0, 1)) + 
    scale_color_discrete( labels = c("Control", "TF2Exp", "Add TF-TF interaction", 'Add HiC'), guide = guide_legend(override.aes = list(size=2), title = NULL))
 
random_tf_plot
hic_plot <- f_plot_performance_and_stats_test(sample_performance_merge, 'TFaddPenalty', 'TF')
plot(hic_plot)
ggsave(f_p('%s/hic_plot.tiff', figure_dir), width =7, height =7, plot = hic_plot)

SNP_All_plot = f_plot_performance_and_stats_test(good_performance, 'SNP', 'All')
f_plot_performance_and_stats_test(good_performance, 'SNP', 'TF', test_stats = T)
f_plot_performance_and_stats_test(good_performance, 'SNP', 'SNPinTF')
f_plot_performance_and_stats_test(good_performance, 'SNP', 'TF.snp')




#f_plot_performance_and_stats_test(sample_performance_merge, 'TF.snp', 'TFsnpMatch')

#f_plot_performance_and_stats_test(sample_performance_merge, 'TF.snp', 'SNPinTF')


f_plot_performance_and_stats_test(sample_performance_merge, 'TF' ,'TFaddInteract', publication = F)
TF_interaction_plot <- f_plot_performance_and_stats_test(sample_performance_merge, 'TF' ,'TFaddInteract', publication = F)
ggsave(filename = f_p('%s/tf_interaction_plot.tiff', figure_dir), plot = TF_interaction_plot, width = 7, height = 7, units = 'in')

compare_rare <- f_plot_performance_and_stats_test(good_performance, 'TFsnpMatch', 'TF.snp' , publication = T)
compare_rare_pub <- compare_rare + xlab('R2 after adding uncommon variants in common-variatnts peaks') +
    ylab('R2 of TF2Exp on common variants')
compare_rare_pub

compare_rare_new_peaks <- f_plot_performance_and_stats_test(good_performance, 'TF', 'TFsnpMatch', publication = T)
compare_rare_new_peaks_pub <- compare_rare_new_peaks + ylab('R2 befroe adding uncommon-variatnts peaks') +
    xlab('R2 after adding uncommon-variant only peaks')
compare_rare_new_peaks


TF_SNPinTF <-f_plot_performance_and_stats_test(good_performance, 'TF', 'TF.snp', publication = T) + xlab('R2 of TF2Exp on all variants') +
    ylab('R2 of TF2Exp on SNPs')



SNP_SNPinTF <-f_plot_performance_and_stats_test(good_performance, 'SNP', 'SNPinTF', publication = T) + xlab('R2 of SNP-based models on all SNPs') +
    ylab('R2 of SNP-based models on all SNPs in TF binding peaks')
SNP_SNPinTF


#Should use good genes, as lower performance genes are meaning less.
#f_plot_performance_and_stats_test(sample_performance_merge, 'SNPinTF', 'TF.snp', publication = F)
#f_plot_performance_and_stats_test(sample_performance_merge, 'SNP', 'TF.snp', publication = F)


f_plot_performance_and_stats_test(good_performance, 'SNPinTF', 'TF.snp', publication = F)
f_plot_performance_and_stats_test(good_performance, 'SNP', 'TF.snp', publication = F)




SNPinTF_TFsnp_plot = f_plot_performance_and_stats_test(good_performance, 'SNPinTF', 'TF.snp', publication = 'Normal')

snp_tf_cmp <- SNPinTF_TFsnp_plot + xlab('Model performance of SNP-based models on SNPinTF data set') + ylab('Model performace of TF2Exp on SNPinTF data set')
snp_tf_cmp


f_plot_performance_and_stats_test(sample_performance_merge, 'SNP', 'random')
plots_snp_tf=f_compare_improvment_for_two_groups('SNP', "TF", good_performance, debug = FALSE)
plots=f_compare_improvment_for_two_groups('SNP', "All", good_performance, debug = FALSE)
plots=f_compare_improvment_for_two_groups('SNPinTF', "TFsnpMatch", good_performance, debug = FALSE)
plots=f_compare_improvment_for_two_groups('SNPinTF', "TF.snp", good_performance, debug = FALSE)

snp_tfsnp_plots = f_compare_improvment_for_two_groups('SNP', "TF.snp", good_performance, debug = FALSE)



TF.SNP_plots=f_compare_improvment_for_two_groups('TF.snp', "TFsnpMatch", good_performance, debug = FALSE)

TF_TF.SNP_plots=f_compare_improvment_for_two_groups("TFsnpMatch", 'TF' ,good_performance, debug = FALSE)

TF_TF.SNP_plots$density


#output:

hic_tf_cmp_stats = f_plot_performance_and_stats_test(sample_performance_merge, 'TFaddPenalty', 'TF', test_stats = T)

cat('Compare addPenalty with TF model in all genes (addPenalty - TF model), p-value (two sided)', hic_tf_cmp_stats$pvalue, 'Median performance diff',
    hic_tf_cmp_stats$median_difference, 'Mean performance diff', hic_tf_cmp_stats$mean_diff, '\n', file = perf_file, append = T)


interaction_tf_cmp_stats = f_plot_performance_and_stats_test(sample_performance_merge, 'TF', 'TFaddInteract', test_stats = T)
cat('Compare TF and TFinteract model in all genes (first - second), p-value (two sided)', interaction_tf_cmp_stats$pvalue, 'Median performance diff',
    interaction_tf_cmp_stats$median_difference, 'Mean performance diff', interaction_tf_cmp_stats$mean_diff, '\n', file = perf_file, append = T)


SNP_TFsnp_stats = f_plot_performance_and_stats_test(good_performance, 'SNP', 'TF.snp', test_stats = T)
cat('Compare SNP and TF.SNP model in well predicted genes (first - second), p-value (two-sided)', SNP_TFsnp_stats$pvalue, 'Median performance diff',
    SNP_TFsnp_stats$median_difference, 'Mean performance diff', SNP_TFsnp_stats$mean_diff, 'Percentage improvment', SNP_TFsnp_stats$mean_diff/SNP_TFsnp_stats$mean2, '\n', file = perf_file, append = T)



SNP_TFsnp_stats = f_plot_performance_and_stats_test(sample_performance_merge, 'SNP', 'TF.snp', test_stats = T)
cat('Compare SNP and TF.SNP model in all genes (first - second), p-value (two-sided)', SNP_TFsnp_stats$pvalue, 'Median performance diff',
    SNP_TFsnp_stats$median_difference, 'Mean performance diff', SNP_TFsnp_stats$mean_diff, '\n', file = perf_file, append = T)


SNPinTF_TFsnp_stats = f_plot_performance_and_stats_test(sample_performance_merge, 'SNPinTF', 'TF.snp', test_stats = T)
cat('Compare SNPinTF with TF.SNP in all genes, p-value', SNPinTF_TFsnp_stats$pvalue, 'Median performance diff',
    SNPinTF_TFsnp_stats$median_difference, 'Mean performance diff', SNPinTF_TFsnp_stats$mean_diff, '\n', file = perf_file, append = T)

SNPinTF_TFsnp_stats = f_plot_performance_and_stats_test(good_performance, 'SNPinTF', 'TF.snp', test_stats = T)
cat('Compare SNPinTF with TF.SNP in shared predictable genes, p-value', SNPinTF_TFsnp_stats$pvalue, 'Median performance diff',
    SNPinTF_TFsnp_stats$median_difference, 'Mean performance diff', SNPinTF_TFsnp_stats$mean_diff, '\n', file = perf_file, append = T)


ggsave(filename = f_p('%s/perf_tf_and_random.tiff', figure_dir), plot = random_tf_plot, width = 7, height = 7, units = 'in')

ggsave(filename = f_p('%s/cmp_tf_snp.tiff', figure_dir), plot = snp_tf_cmp, width = 7, height = 7, units = 'in')

ggsave(filename = f_p('%s/rare_new_peak_plot.tiff', figure_dir), plot = compare_rare_new_peaks_pub, width = 7, height = 7, units = 'in')
ggsave(filename = f_p('%s/rare_plot.tiff', figure_dir), plot = compare_rare_pub, width = 7, height = 7, units = 'in')

rare_plots = arrangeGrob( compare_rare_pub, compare_rare_new_peaks_pub, ncol =2)
ggsave(filename = f_p('%s/rare_plots.tiff', figure_dir), plot = rare_plots, width = 14, height = 7, units = 'in')



###########################Pairwise comparisn done#################################



################Section 5: Test the gene mean and var against the performance###########
batch_expression = read.table(f_p('./data/%s/rnaseq/GEUVADIS.Gene.DATA_MATRIX', loc_batch_name), header = TRUE)
rownames(batch_expression) = batch_expression$gene
batch_expression$gene = NULL
#batch_expression_log = log(batch_expression + 1)
batch_expression_log = batch_expression
expression_array=data.frame( gene = rownames(batch_expression),  mean = rowMeans(batch_expression_log, na.rm = T), var = apply(batch_expression_log, 1, var, na.rm = T))

rownames(expression_array) = expression_array$gene
head(expression_array)

#head10(expression_array)
tf_performance_merge <- sample_performance_merge %>% filter(mode == 'TF')
head(tf_performance_merge)
tf_performance_merge$mean_exp = expression_array[tf_performance_merge$gene,'mean']
tf_performance_merge$var = expression_array[tf_performance_merge$gene,'var']
plot_performance = subset(tf_performance_merge, performance > threshold)

f_plot_values_against_performance(plot_performance, 'mean_exp')
perf_var_plot <- f_plot_values_against_performance(plot_performance, 'var')
perf_var_stats <- f_plot_values_against_performance(plot_performance, 'var', test_stats =T)
perf_var_stats

cat('Corelate the performance of gene and variance of gene expression, pvalue:', perf_var_stats$pvalue, 'Cor-coef:',
    perf_var_stats$coef, '\n', file = perf_file, append = T)

var_plot <- perf_var_plot + theme_Publication(12) + xlab('Model performance') + ylab('Variance of gene expression')
plot(var_plot)
ggsave(filename = f_p('%s/var_plot.tiff', figure_dir), plot = var_plot, width = 7, height = 7, units = 'in')
#####################gene Var done######################################




###############Section 6: GO analysis, wait until the whole genome analysis.##################
source('s_GO_analysis_predictions.R')
#thres is working.
f_performance_go_analysis(sample_performance_merge, 'TF', thres = 0.59) #check this after the whole genome analysis.
#########################################
all_gene_loc = read.table('./data/raw_data/rnaseq/all.ensemble.genes.gene_start', header = T)

head(all_gene_loc)

rownames(all_gene_loc) = all_gene_loc$ensembl_gene_id

great_genes <- sample_performance_merge %>% filter(mode == 'TF', performance > 0.4)%>%
    arrange(-performance) %>% dplyr::select(gene, performance)
head(great_genes)
dim(great_genes)
head(all_gene_loc)
great_bed = all_gene_loc[str_replace(great_genes$gene, '[.].*', ''), c(1:5,7) ]
head(great_bed)
colnames(great_bed) = c('gene', 'chr', 'start', 'end', 'strand', 'name')
great_bed$tss_start = great_bed$start
great_bed$tss_start[ great_bed$strand == '-1' ] = great_bed[ great_bed$strand == '-1', 'end' ]

great_bed$chr = paste0('chr', great_bed$chr)
great_bed$name = paste0(great_bed$gene, ':', great_bed$name)
head(great_bed)

bed_data <- great_bed %>% mutate( bed_start = as.numeric(tss_start) - 500, bed_end = as.numeric(tss_start) + 500 ) %>% dplyr::select(chr, bed_start, bed_end)
write.table(great_bed[1:400,c( 'chr', 'start', 'end', 'gene')], file = f_p('%s/genes_for_great.bed', results_dir), quote = F, sep = '\t', row.names = F, col.names = F)
write.table(great_bed[1:400,c( 'chr', 'start', 'end', 'name')], file = f_p('%s/genes_for_great2.bed', results_dir), quote = F, sep = '\t', row.names = F, col.names = F)


stop()

################The Q-Q plot#######
if (FALSE){
    f_plot_performance_and_stats_test(good_performance, 'All', 'TF')
    f_plot_performance_and_stats_test(good_performance, 'SNPinTF', 'TF')
    f_plot_performance_and_stats_test(good_performance, 'TFsnpMatch', 'TF')
    f_plot_performance_and_stats_test(good_performance, 'TFaddInteract', 'TF')
    f_plot_performance_and_stats_test(good_performance, 'SNPinTF', 'SNP')
    f_compare_improvment_for_two_groups('SNP', "All", good_performance, thres = 0.01)
    f_compare_improvment_for_two_groups('SNP', "TF", good_performance)
    f_compare_improvment_for_two_groups('SNP', "SNPinTF", good_performance, thres = 0.01)
    f_compare_improvment_for_two_groups('SNPinTF', "SNP", good_performance, thres = 0.01)

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
}




#' #Test performance drop against model variance
try(f_test_performance_drop_with_ml_variance(sample_performance_merge, 'SNP', 'SNPinTF'))

sample_performance_merge$Model = str_replace(sample_performance_merge$mode, 'random', 'Random')

ggplot(subset(sample_performance_merge, gene %in% intersect_genes & mode %in% c('TF', 'random' ,'Random')), aes(x = performance , color = mode )) + geom_density()




#' Analysis the performance in the external datasets.
head(sample_performance_merge)

plot_performance <- good_performance %>% filter( mode %in% c('TF', 'SNP'))

table(plot_performance$mode)

library(tidyr)
external_cmp <- plot_performance %>% select(gene, performance, mode) %>% spread(value = performance, key = mode )


f_debug <- function(){
    source('s_summary_fun.R')
    sample_performance_merge = f_collect_performance_in_multiple_chrs(loc_batch_name, c(mode_list[1]), chr_num_list)
}




