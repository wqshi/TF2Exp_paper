#This script is for the feature analysis part.
source('~/R/s_ggplot2_theme.R')
source('s_summary_fun.R')
source('s_project_funcs.R')
source('s_gene_regression_fun.R')
source('r_feature_analysis_fun.R')
library(optparse)
library(dplyr)
library(xtable)
option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="dataset file name", metavar="character"),
    make_option(c("--collection"),      type="character", default='addPenalty', help="The parameter based name", metavar="character"),
    make_option(c("--chr_batch"),      type="character", default='1chr', help="1chr for chr22, 3chr for (2, 10, 22)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
str(opt)

if (is.null(opt$batch_name)){
    #batch_name = '54samples_evalue'
    #batch_name = '462samples_sailfish_quantile'
    #batch_name = '462samples_snyder_norm'
    #batch_name = '462samples_quantile'
    #batch_name = '462samples_quantile_rmNA'
    #batch_name = '462samples_log_quantile'
    loc_batch_name = '358samples_regionkeepLow'
    ##chr_str = 'chr22'
    chr_num_list = c(10, 15, 22)
    collection_name = 'peer358'
    #collection_name = 'addPenalty'
}else{
    loc_batch_name = opt$batch_name
    if (opt$chr == '1chr'){
        chr_num_list = c(22)
    }else{
        chr_num_list = c(22, 10, 2)
    }
    collection_name = opt$collection
}

threshold = 0.05

#####Count the feature categories######
mode_list = modes_list[[collection_name]]
mode_list = f_create_new_mode_list(mode_list, 'trad', 'cor')

mode_list


#' #Check the sign of the coefficient in promoter and enhancer regions
features_stats = f_merge_selected_features(loc_batch_name, paste0('chr', chr_num_list), mode_list['TF'], debug = F)

features_merge = subset(features_stats$features, !is.na(feature_start))

features_merge$tf_rename = str_replace(features_merge$rename, '(enhancer|promoter)', '')
table(features_merge$tf_rename)





well_genes = features_stats$performances %>% filter(performance > 0.05) %>% select(gene)
features_merge$location = 'promoter'
features_merge$location[ grepl('enhancer', features_merge$rename) ] = 'enhancer'

well_features_merge = features_merge %>% filter(gene %in% well_genes$gene)
table(well_features_merge$location)

head(well_features_merge)

promoter_coef = subset(well_features_merge, location == 'promoter')$score
enhancer_coef = subset(well_features_merge, location == 'enhancer')$score
wilcox_test = wilcox.test( promoter_coef , enhancer_coef)



##Output: Final plot to compare the coeficients in promoters and enhancers.
coef_box_plot<- ggplot(well_features_merge, aes(location, score )) + geom_boxplot(notch = TRUE) +
    ylab('Feature coefficients') + xlab('Feature location')+ annotate("text", label = f_p('p-value: %.1e\n enhacer mean %.2f\n promoter mean %.2f ',
                                                                   wilcox_test$p.value, mean(enhancer_coef), mean(promoter_coef) ),
                 x = 1.2, y = 0.4, size = 6, colour = "red") + theme_Publication(base_size = 12)


ggplot(well_features_merge, aes(location, score )) + geom_boxplot(notch = TRUE) +
    ggtitle('Feature coefficients') + annotate("text", label = f_p('p-value: %.1e\n enhacer mean %.2f\n promoter mean %.2f ',
                                                                   wilcox_test$p.value, mean(enhancer_coef), mean(promoter_coef) ),
                 x = 1.2, y = 0.4, size = 6, colour = "red")


ggplot(well_features_merge, aes(score, ..density.., color = location )) + geom_freqpoly() +
    ggtitle('Feature coefficients') + annotate("text", label = f_p('p-value: %.1e\n enhacer mean %.2f\n promoter mean %.2f ',
                                                                   wilcox_test$p.value, mean(enhancer_coef), mean(promoter_coef) ),
                 x = -0.4, y = 12, size = 6, colour = "red")










############ Functional role of TF binding towards gene regulation#####################
feature_summary_stats <- features_merge %>% group_by(rename) %>% dplyr::summarise(positive = sum(score > 0), total = length(score), pos_ratio = positive/total, mean_score = mean(score)) %>%
    arrange(desc(pos_ratio)) %>% filter(total > 8) %>% arrange(pos_ratio)
head(well_features_merge)

#tf_stats = well_features_merge %>% group_by(rename) %>% dplyr::summarise(mean_score = mean(score))
tf_stats = feature_summary_stats
tf_stats$location = 'promoter'
tf_stats$location[grepl('enhancer',tf_stats$rename)] = 'enhancer'
tf_stats$tf = str_replace(tf_stats$rename, 'promoter|enhancer', '')



#' Output: Top 3 negative and Top 3 positive
head(tf_stats, n = 5)
tail(tf_stats, n = 5)




#colnames(collect_table_write_sort) = c('Feature', 'Location', "Selected",  "Unselected",  "P-value", "Odds ratio")
#print.xtable(xtable(collect_table_write_sort), type = 'html', './result/figures/feature.html', include.rownames=FALSE)


#Number of TFs show positive coefficient.
feature_summary_stats %>% summarise( n = sum(pos_ratio > 0.5)/length(pos_ratio) )

write.table(format(tf_stats, digits = 3), file = './result/stats/pos_ratio.txt', quote = F, row.names = F)



library(tidyr)

head(tf_stats)

tf_stats$pos_ratio

plot_data <- tf_stats %>% select(tf, location, pos_ratio) %>% spread(key = location, value = pos_ratio)

head(plot_data)
pos_ratio_plot <- ggplot(plot_data, aes(promoter, enhancer)) + geom_point() + geom_abline() +
    geom_text(data = plot_data[complete.cases(plot_data),], aes(x = promoter, y = enhancer + 0.02, label = tf)) +
    xlab('Positive ratio in promoters') + ylab('Positive ratio in enhancer') + theme_Publication(12)

plot(pos_ratio_plot)



#' #Calculate the mean distance of the features in promter and enhancer regions.
tf_bed = read.table('./data/358samples_regionkeepLow/rnaseq/transcript_data.bed', header = T)
gene_bed = tf_bed[,1:6]

gene_bed$tss = gene_bed$start + 2000
rownames(gene_bed) = gene_bed$gene
head(gene_bed)

features_merge$tss = gene_bed[features_merge$gene,'tss']
features_merge$feature_dis = features_merge$feature_start - features_merge$tss

head(features_merge)

ggplot(features_merge, aes(log10(abs(feature_dis)), color = location)) + geom_density()

features_merge$upstream = 'Upstream'
features_merge$upstream[features_merge$feature_dis > 0] = 'Downstream'

##Todo: Add the overall background
ggplot(features_merge, aes(log10(abs(feature_dis)), color = location)) + geom_density() + facet_wrap(~upstream, nrow =2)

features_merge %>% group_by(location) %>% dplyr::summarise(mean_dis = mean(abs(feature_dis)))

features_merge$positive = 'positive'
features_merge$positive[features_merge$score < 0] = 'negative'
ggplot(features_merge, aes(abs(score), abs(feature_dis))) + geom_point()


##Output
#distal_plot <- ggplot(features_merge, aes(sign(feature_dis)*log10(abs(feature_dis)), color = location)) + geom_freqpoly() + theme_Publication(12)

ggplot(features_merge, aes(feature_dis)) + stat_ecdf(geom = "step") + theme_Publication(12) + facet_wrap(~location, nrow =2, scales = 'free_x')
distal_plot <- ggplot(features_merge, aes(feature_dis)) + geom_histogram() + theme_Publication(12) + facet_wrap(~location, nrow =2, scales = 'free_x')

library(gridExtra)
feature_coef_plot <- arrangeGrob( pos_ratio_plot, distal_plot, ncol = 2)

ggsave('./result/figures/feature_coef.tiff', plot = feature_coef_plot, width = 14, height = 7)

positive_location_test=fisher.test( as.matrix(table(features_merge$positive, features_merge$upstream)))



#Get the top genes, and get the features.
performance_df = features_stats$performances
performance_df_sorted = f_sort_by_col(performance_df, 'performance', TRUE)
head(performance_df_sorted)

feature_df = features_stats$features
head(feature_df)
top_gene_features = feature_df %>% filter(gene == performance_df_sorted$gene[1])

top_gene_features


stop()

#################Section: Stats of the selected features.################
collect_table = f_feature_enrichment( loc_batch_name, c('chr22', 'chr10', 'chr15'), mode_list['TF'], debug = TRUE)
#f_feature_enrichment( loc_batch_name, c('chr22'), mode_list['All'], debug = TRUE)
#f_feature_enrichment( loc_batch_name, c('chr22', 'chr2', 'chr10'), mode_list['TF'])

head(collect_table)
colnames(collect_table) = c('feature', 'selected', 'unselected', 'p-value', 'Odds_ratio')

collect_table$location = 'Enhancer'
collect_table$location[grep('promoter', collect_table$feature,)]  = 'Promoter'
collect_table$location[grep('P_E', collect_table$feature,)]  = 'P_E-Interaction'
collect_table$location[grep('TFoverlap', collect_table$feature,)]  = 'TF-interaction'
collect_table$feature = str_replace(collect_table$feature, '(promoter|enhancer|P_E|TFoverlap[.])', '')
table(collect_table$location)
collect_table_write = collect_table
collect_table_write['p-value'] = format(as.numeric(collect_table[,'p-value']), digits = 2)
collect_table_write['Odds_ratio'] = format(as.numeric(collect_table[,'Odds_ratio']), digits = 2)
head(collect_table_write)
colnames(collect_table_write)

collect_table_write = collect_table_write[,c('feature', 'location', "selected",  "unselected",  "p-value", "Odds_ratio")]

collect_table_write$selected = as.numeric(collect_table_write$selected)

collect_table_write_sort = f_sort_by_col(collect_table_write, 'selected', TRUE)
#install.packages('ReporteRs')
#install.packages('gdtools')
#f_table_to_word(collect_table_write)

collect_table_write_sort

View(subset(collect_table_write_sort))
#install.packages('xtable')
library('xtable')

colnames(collect_table_write_sort) = c('Feature', 'Location', "Selected",  "Unselected",  "P-value", "Odds ratio")
print.xtable(xtable(collect_table_write_sort), type = 'html', './result/figures/feature.html', include.rownames=FALSE)

collect_table_write %>% filter(location == 'Promoter')


##Correlation between number of peaks and selected features.
peak_count = read.table('~/expression_var/R/data/raw_data/tf/encode_peaks/processed/peak.counts')
colnames(peak_count) = c('count', 'file_name')

peak_count$tf = toupper(str_replace(str_replace(peak_count$file_name,'.*gm12878-', ''), '.narrowPeak',''))
rownames(peak_count)=peak_count$tf

f_feature_count_correlation <- function(input_table, peak_count){
    cor.test(input_table$selected,   peak_count[input_table$feature,'count'])
}

#Features peak count and selected number of features.
f_feature_count_correlation(subset(collect_table_write, location == 'Promoter'), peak_count)
f_feature_count_correlation(subset(collect_table_write, location == 'Enhancer'), peak_count)


#The median odds ratio:
sum(collect_table_write_sort[1:10,'Selected'])/sum(collect_table_write_sort[,'Selected'])
collect_table %>% group_by(location) %>% dplyr::summarise(mean_odds = median(as.numeric(Odds_ratio)))


if (FALSE){

chr_list=paste0('chr', chr_num_list)


#Count the number of features in the final selected features.
all_feature_table = f_summary_selected_features(loc_batch_name, chr_list, mode_list['All'], output_feature_file = 'feature.table', threshold, return_tf = T)
SNP_feature_table = f_summary_selected_features(loc_batch_name, c('chr22'), mode_list['SNP'], output_feature_file = 'feature.table.tf', threshold )
tf_feature_table = f_summary_selected_features(loc_batch_name, chr_list, mode_list['TF'], output_feature_file = 'feature.table.tf', threshold, return_tf = T)
tf_feature_table = f_summary_selected_features(loc_batch_name, chr_list, mode_list['TF'], output_feature_file = 'feature.table.tf.location', threshold, return_tf = F)
head(all_feature_table)
nrow(all_feature_table)




##Compare the difference between All and SNP model. The number of selected features and lambda.
return_list1 = f_summary_regression_results(loc_batch_name, 'chr22', mode_list['All'], rsync_flag = FALSE, return_features = TRUE)
return_list2 = f_summary_regression_results(loc_batch_name, 'chr22', mode_list['SNP'], rsync_flag = FALSE, return_features = TRUE)

head(return_list1$features)

library('stringr')
gene_list = unique(str_replace(return_list1$features$name, '[|].*', ''))
gene_list2 = unique(str_replace(return_list2$features$name, '[|].*', ''))

shared_list = intersect(gene_list, gene_list2)
loc_gene = 'ENSG00000130489.8'
collect_stats = data.frame()
for (loc_gene in shared_list){
    if( return_list1$performance[loc_gene, 'performance'] > 0.05 & return_list2$performance[loc_gene, 'performance'] > 0.05 ){
       
        feature_set1=return_list1$features[grep(loc_gene, return_list1$features$name),]
        feature_set2=return_list2$features[grep(loc_gene, return_list2$features$name),]
        
        if( all(unique(feature_set1$rename) %in% c('SNP', '(Intercept)')) & all(unique(feature_set2$rename) %in% c('SNP', '(Intercept)')) ){
            cat(loc_gene,'All', return_list1$performance[loc_gene, 'performance'],
                 'SNP', return_list2$performance[loc_gene, 'performance'],
                'Set diff', length(setdiff(feature_set1$name, feature_set2$name)),
                nrow(feature_set1), nrow(feature_set2),
                return_list1$performance[loc_gene, 'lambda'],
               return_list2$performance[loc_gene, 'lambda'],'\n')
            collect_stats = rbind(collect_stats, c(return_list1$performance[loc_gene, 'performance'], return_list2$performance[loc_gene, 'performance']))
        }
    }
}

colnames(collect_stats) = c('All', 'SNP')
collect_stats =as.data.frame( data.matrix(collect_stats))
apply(collect_stats, 2, mean)
head(collect_stats)
class(collect_stats)
sum(collect_stats$All < collect_stats$SNP)
sum(collect_stats$All == collect_stats$SNP)
nrow(collect_stats)

    
#######Find top 50 genes####################
all_performance <- sample_performance_merge %>% filter(mode == 'All', chr == 'chr22')
best_gene_list=(f_sort_by_col(all_performance, 'performance', decr_flag = T)[1:50,])$gene



################Best genes#####################
best_index=which.max(sample_performance_merge$performance)
best_gene = sample_performance_merge[best_index,'gene']
sample_performance_merge %>% filter(gene == best_gene)


feature_table = return_list2$features
feature_table[grepl(best_gene, feature_table$name),]


##############Sort features#################
    for (best_gene in best_gene_list){

        all_best_gene = f_get_feature_file_for_one_gene(best_gene, mode_list['All'])
        SNP_best_gene = f_get_feature_file_for_one_gene(best_gene, mode_list['SNP'])

        all_features = rownames(all_best_gene)
        SNP_features = rownames(SNP_best_gene)
        shared_cols=intersect(all_features, SNP_features)
        all_unique = setdiff(all_features, SNP_features)
        SNP_unique = setdiff(SNP_features, all_features)

        print(f_p('shared %s, all_unique %s, SNP_unique %s', length(shared_cols), length(all_unique), length(SNP_unique)))
                                        #print(SNP_best_gene[SNP_unique,])
                                        #print(all_best_gene[all_unique,])
    }

    f_debug <- function(){

        source('s_summary_fun.R')
        all_feature_table = f_summary_selected_features(loc_batch_name, c('chr22'), mode_list['All'], output_feature_file = 'feature.table', threshold )

    }
}



#Step###########Overlap with the FANTOM5 enhancers#############
library(dplyr)
library(tidyr)
source('r_bedtools.R')


##Read f5 enhancer-promoter, transfer it to enhancer-gene_id pairs.
fantom5_cor = read.table('./data/raw_data/fantom5/hg19_enhancer_promoter_correlations_distances_cell_type.txt', header = T)
f5_promoter <- fantom5_cor %>% select(promoter, enhancer) %>% separate(promoter, into = c('chr', 'start', 'end', 'strand'))
tss_df = read.table('./data/raw_data/rnaseq/transcript_loc.bed', header = T, sep = '\t')

overlapped_bed = f_bedtools(f5_promoter[,c('chr', 'start', 'end', 'enhancer')], tss_df[,c('chr', 'start', 'end', 'gene')], fun = 'intersect', paras = '-wao')
enhancer_gene <- overlapped_bed %>% dplyr::rename( enhancer = V4, gene = V8) %>% dplyr::select(enhancer, gene) %>%
     tidyr::separate(enhancer, into = c('chr', 'start', 'end')) %>% distinct() %>% arrange(gene)

enhancer_gene_filter = enhancer_gene %>% filter(gene != '.')



##Overlap selected features of each gene with the enhancer-gene pairs/ compare with unselected features.
f_overlap_tf_features_with_f5_enhancers <- function(feature_bed, enhancer_gene_filter, good_gene_list){
    feature_enhancer_bed = feature_bed[grepl('enhancer', feature_bed$name), ]
    dim(feature_enhancer_bed)
    head(enhancer_gene_filter)

    overlap_df = f_bedtools(feature_enhancer_bed, enhancer_gene_filter, paras = '-wao')
    head(overlap_df)
    overlap_df_sub = as.data.frame(overlap_df)[c('V1', 'V2', 'V3', 'V4', 'V8')]

    head(overlap_df_sub)
    colnames(overlap_df_sub) =c('chr','start', 'end', 'name', 'matched_gene' )
    overlap_df_gene  = overlap_df_sub %>% mutate( gene = str_replace(name, '[.].*', '') ) %>% filter( matched_gene == gene | matched_gene == '.', gene %in% good_gene_list)

    print(table(overlap_df_gene$gene == overlap_df_gene$matched_gene))
    return (overlap_df_gene)
}

unselected_features = features_stats$merge_control
selected_features = features_stats$features


good_performance_gene <- features_stats$performances %>% filter(performance > 0.05)

selected_bed = f_extract_bed_from_features(selected_features) 
control_bed = f_extract_bed_from_features(unselected_features)

good_gene_list = str_replace(good_performance_gene$gene, '[.][0-9]*', '')
real_overlap = f_overlap_tf_features_with_f5_enhancers(selected_bed, enhancer_gene_filter, good_gene_list)
control_overlap = f_overlap_tf_features_with_f5_enhancers(control_bed, enhancer_gene_filter, good_gene_list)v

