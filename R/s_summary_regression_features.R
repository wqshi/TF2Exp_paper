
#This script is for the feature analysis part.
source('~/R/s_ggplot2_theme.R')
source('s_summary_fun.R')
source('s_project_funcs.R')
source('s_gene_regression_fun.R')
source('r_feature_analysis_fun.R')
library(optparse)
library(dplyr)
library(xtable)
library(tidyr)
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
    chr_num_list = 1:22
    collection_name = 'peer358corRmdup'
    #collection_name = 'elastic'
    #collection_name = 'peer358cor'
    #collection_name = 'rmZero'
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
source('./s_test_cat.R')



#######################1 Prepare the feature data###############################
mode_list = modes_list[[collection_name]][c('TF', 'SNP')]
mode_list = f_create_new_mode_list(mode_list, 'trad', 'cor')
read_flag = FALSE
results_dir = f_p('./data/%s/rnaseq/results/%s/', loc_batch_name, collection_name)
figure_dir = f_p('%s/figures/', results_dir)
dir.create(figure_dir)
output_file = f_p('%s/feature_results.txt', results_dir)
f_write_paper_results('=====Features numbers======', data = date(), file = output_file)


if (read_flag == TRUE){
    features_stats = f_merge_selected_features(loc_batch_name, paste0('chr', chr_num_list), mode_list['TF'], debug = F)
    save(features_stats, file = f_p('%s/features_stats', results_dir))
}else{
    load(file = f_p('%s/features_stats', results_dir))
}

features_merge = subset(features_stats$features, !is.na(feature_start))
features_merge$tf_rename = str_replace(features_merge$rename, '(enhancer|promoter)', '')
table(features_merge$tf_rename)
features_merge$location = 'promoter'
features_merge$location[ grepl('enhancer', features_merge$rename) ] = 'enhancer'

features_merge$tf_rename[features_merge$tf_rename == 'PU'] = 'PU.1'
features_merge$tf_rename[features_merge$tf_rename == 'USF'] = 'USF.1'

sample_performance_merge = read.table(file = f_p('%s/sample_performance_merge', results_dir), header = T)



####Write key features and performance of the gene######
features_write = features_stats$features
features_write$tf_rename = str_replace(features_write$rename, '(enhancer|promoter)', '')
features_write$location = 'promoter'
features_write$location[ grepl('enhancer', features_write$rename) ] = 'distal'
features_write$location[ is.na(features_write$feature_start)] = NA

features_write$tf_rename[features_write$tf_rename == 'PU'] = 'PU.1'
features_write$tf_rename[features_write$tf_rename == 'USF'] = 'USF.1'

features_write_select = features_write %>% select(gene, chr, feature_start, feature_end, tf_rename, location, score) %>%
    rename(tf = tf_rename)

performance_write = features_stats$performances %>% select(gene, train_performance, performance, alpha, lambda)
head(performance_write)
git_data_dir = '../results/'
write.table(features_write_select, file = f_p('%s/features_all', git_data_dir),quote = FALSE, sep='\t', row.names = FALSE)
write.table(performance_write, file = f_p('%s/performance_all', git_data_dir),quote = FALSE, sep='\t', row.names = FALSE)

head(features_stats$performances)



#################: Compare the coefficient in promoter and enhancer regions#############
str(features_stats, max.level = 1)

well_genes = features_stats$performances %>% filter(performance > 0.05) %>% dplyr::select(gene)

well_features_merge = features_merge %>% filter(gene %in% well_genes$gene, abs(score) > 1e-6)
table(well_features_merge$location)

head(well_features_merge)

well_features_merge %>% group_by(gene, tf_rename, feature_start, feature_end) %>% nrow()
hist(well_features_merge$score, breaks= 100)

well_features_merge %>% filter(abs(score) < 1e-4) %>% nrow()

#Test the coef difference in promoter and distal regions.
promoter_coef = subset(well_features_merge, location == 'promoter')$score
enhancer_coef = subset(well_features_merge, location == 'enhancer')$score
wilcox_test = wilcox.test( abs(promoter_coef) , abs(enhancer_coef), alternative = 'two.sided', conf.int = T)

f_write_paper_results('Test abs(promoter_coef) > abs(enhancer_coef), p-value:', wilcox_test$p.value, output_file)
f_write_paper_results('Mean promoter coef:', mean(abs(promoter_coef)), output_file)
f_write_paper_results('Mean enhancer coef:', mean(abs(enhancer_coef)), output_file)
f_write_paper_results('Estimated median enhancer coef:', wilcox_test$estimate, output_file)
f_write_paper_results('Number of loction:', table(well_features_merge$location), output_file)


#Final plot to compare the coeficients in promoters and enhancers.
well_features_merge$location2 = 'Distal'
well_features_merge[well_features_merge$location == 'promoter', 'location2'] = 'Promoter'
head(well_features_merge)


coef_box_plot<- ggplot(well_features_merge, aes(location2, abs(score) )) + geom_boxplot() +
    ylab('Absolute feature effect sizes') + xlab('Feature location')+ annotate("text", label = f_p('P-value: %.1e',
                                                                   wilcox_test$p.value),
                 x = 1.5, y = 0.7) + theme_Publication(base_size = 12) +
    geom_hline(yintercept=0, linetype = '1F')

plot(coef_box_plot)
ggsave(f_p('%s/feature_coef_boxplot.tiff', figure_dir), plot = coef_box_plot)

hist(well_features_merge$score, breaks = 100)







############### Functional role of TF binding towards gene regulation#####################
######Positive ratio in proximal and distal regions for each TF###########

##Prepare the tf-based data
tf_stats <- features_merge %>% group_by(rename) %>% dplyr::summarise(positive = sum(score > 0), total = length(score), pos_ratio = positive/total, mean_score = mean(score)) %>%
    arrange(desc(pos_ratio)) %>% filter(total > 10) %>% arrange(pos_ratio)
tf_stats$location = 'promoter'
tf_stats$location[grepl('enhancer',tf_stats$rename)] = 'enhancer'
tf_stats$tf = str_replace(tf_stats$rename, 'promoter|enhancer', '')


## Top 3 negative and Top 3 positive
f_write_paper_results('Top 5 negative TFs', tf_stats[1:5,], output_file )
f_write_paper_results('Top 5 positive TFs', tf_stats[(nrow(tf_stats) -5):nrow(tf_stats),], output_file )
head(tf_stats, n = 5)
tail(tf_stats, n = 5)


#Number of TFs show positive coefficient.
tf_stats %>% summarise( n = sum(pos_ratio > 0.5)/length(pos_ratio) )
write.table(format(tf_stats, digits = 3), file = f_p('./%s/pos_ratio.txt', results_dir), quote = F, row.names = F)


plot_data <- tf_stats %>% dplyr::select(tf, location, pos_ratio) %>% spread(key = location, value = pos_ratio) %>% as.data.frame
label_data<-plot_data[complete.cases(plot_data),]
rownames(label_data) = label_data$tf
extreme_points <-  label_data %>% gather(new, new2, enhancer:promoter) %>% arrange(-new2) %>% filter(row_number() <= 5 | row_number() >= n() - 4)

label_data2 = f_adjust_lable_positions( label_data[extreme_points$tf,] , x_lab = 'promoter', y_lab = 'enhancer', rep.fact = 30,  rep.dist.lmt = 10,
                                      attr.fact = 0.3, iter.max = 20000)
head(label_data2,n =15)
down_tfs = c('SMC3', 'STAT1', 'BCLAF1')
label_data2[down_tfs, 'enhancer'] = label_data2[down_tfs, 'enhancer'] - 0.04
label_data2['BCLAF1', 'promoter'] = label_data2['BCLAF1', 'promoter'] + 0.03
point_color = 'black'
point_color = '#00BFC4'
point_color = 'skyblue4'

pos_ratio_plot <- ggplot(plot_data, aes(promoter, enhancer)) + geom_point(color = point_color) + geom_abline(alpha = 0.3) +
    geom_text(data = label_data2, aes(x = promoter, y = enhancer + 0.02, label = tf)) + xlim(c(0.15, 0.95)) + ylim(c(0.15, 0.95)) +
    xlab('Positive ratio in promoters') + ylab('Positive ratio in distal regulatory regions') + theme_Publication(12)
pos_ratio_plot

#Test: TF in promoter are more positive than distal regions
cor.test(plot_data$promoter, plot_data$enhancer, 'g')
promoter_test = wilcox.test(plot_data$promoter, plot_data$enhancer, 'two.sided')
print(promoter_test)
f_write_paper_results('Test more positive in promoter than in enhancer, p-value:', promoter_test$p.value, output_file)
##Not selected in the manuscript. As the role of TF is complex determined by context.
###################################################################################








#####################Feature distance to gene start positions######################################
#' #Calculate the mean distance of the features in promter and enhancer regions.
tf_bed = read.table('./data/358samples_regionkeepLow/rnaseq/transcript_data.bed', header = T)
gene_bed = tf_bed[,1:6]

gene_bed$tss = gene_bed$start + 2000
rownames(gene_bed) = gene_bed$gene
head(gene_bed)

features_merge$tss = gene_bed[features_merge$gene,'tss']
features_merge$feature_dis = features_merge$feature_start - features_merge$tss

features_merge$upstream = 'Upstream'
features_merge$upstream[features_merge$feature_dis > 0] = 'Downstream'

##Todo: Add the overall background
#ggplot(features_merge, aes(log10(abs(feature_dis)), color = location)) + geom_density() + facet_wrap(~upstream, nrow =2)

features_merge %>% group_by(location) %>% dplyr::summarise(mean_dis = mean(abs(feature_dis)))

features_merge$positive = 'positive'
features_merge$positive[features_merge$score < 0] = 'negative'
head(features_merge)
features_merge$location_rename = 'Promoter'
features_merge$location_rename[features_merge$location == 'enhancer'] = 'Distal regulatory regions'

features_merge$location_rename = factor(features_merge$location_rename, levels = c('Promoter', 'Distal regulatory regions'))

head(features_merge)
table(features_merge$location)

#add fake line to balance the x distribution in promoter
add_fake_line<-features_merge %>% filter(location == 'promoter', feature_dis > 29000)
add_fake_line$feature_dis = -27000
add_fake_line$score = 0


str(features_merge)
features_merge$Score = features_merge$score

features_merge2 = features_merge[,c('Score', 'feature_dis')]
colnames(features_merge2) = colnames(failwith)

str(features_merge2)

colSums(is.na(features_merge2))



distal_plot<-
    ggplot(data=features_merge, aes(y=score, x=feature_dis)) + geom_point(alpha = 0.1, color=point_color, size = 1) +
    geom_density2d(bins = 3, colour = 'green')+
    facet_wrap(~location_rename, nrow = 2, scales = 'free_x') +# scale_size(range = c(0, 0.3)) + 
    theme_Publication(12) + xlab('Feature distance to the gene start site') + ylab('Feature effect size') 

ggsave(f_p('%s/feature_coef_promoter_enhaner.tiff', figure_dir), plot = distal_plot, width =7, height = 7, dpi = 72)
##Output
#distal_plot <- ggplot(features_merge, aes(sign(feature_dis)*log10(abs(feature_dis)), color = location)) + geom_freqpoly() + theme_Publication(12)

positive_location_test=fisher.test( as.matrix(table(features_merge$positive, features_merge$upstream)))

#features_merge_subset <- features_merge %>% filter( ! (feature_dis > 20000 & location_rename == 'Promoter'))
distal_plot2 <-ggplot(features_merge, aes(feature_dis)) + geom_density() + theme_Publication(12) + facet_wrap(~location_rename, nrow =2, scales = 'free') +
    xlab('Feature distance to the gene start site')


library(gridExtra)
tiff(f_p('./%s/feature_coef.tiff', figure_dir), res = 310, width = 10, height = 5, units='in')
#feature_coef_plot <- arrangeGrob(distal_plot, pos_ratio_plot, ncol = 2)
grid.arrange(arrangeGrob(distal_plot, pos_ratio_plot, ncol = 2))
grid.text(label = '(A)',x=unit(0.02, "npc"), y=unit(0.98, "npc"), gp=gpar(fontsize=12))
grid.text(label = '(B)',x=unit(0.52, "npc"), y=unit(0.98, "npc"), gp=gpar(fontsize=12))
dev.off()

tiff(f_p('./%s/feature_distance_coef.tiff', figure_dir), res = 310, width = 10, height = 5, units='in')
#feature_coef_plot <- arrangeGrob(distal_plot, pos_ratio_plot, ncol = 2)
grid.arrange(arrangeGrob(distal_plot, coef_box_plot, ncol = 2))
grid.text(label = '(A)',x=unit(0.02, "npc"), y=unit(0.98, "npc"), gp=gpar(fontsize=12))
grid.text(label = '(B)',x=unit(0.52, "npc"), y=unit(0.98, "npc"), gp=gpar(fontsize=12))
dev.off()



###Output:
ggsave(f_p('%s/feature_coef_promoter_enhaner.tiff', figure_dir), plot = distal_plot, width =7, height = 7)
f_write_paper_results( 'Positive ratio for features in  promoters and distal regulatory regions:',
                      well_features_merge %>% group_by(location) %>% dplyr::summarise( sum(score > 0)/length(score)),
                      output_file)

#######################################################





##############Negative features###########################
features_merge$role = 'positive'
features_merge$role[features_merge$score < 0] = 'negative'

ggplot(subset(features_merge, location == 'promoter' ), aes(abs(feature_dis), color = role)) + geom_density()

#wilcox.test( abs( subset(features_merge, location ==) ) )

##########################################


########Get the top genes, and get the features.################
performance_df = features_stats$performances
performance_df_sorted = f_sort_by_col(performance_df, 'performance', TRUE)
head(performance_df_sorted)

feature_df = features_stats$features
head(feature_df)
top_gene_features = feature_df %>% filter(gene == performance_df_sorted$gene[1])

top_gene_features
######################################################













################Overlap with the FANTOM5 enhancers################################
source('r_bedtools.R')
fantom5_cor = read.table('./data/raw_data/fantom5/hg19_enhancer_promoter_correlations_distances_cell_type.txt', header = T)
f5_promoter <- fantom5_cor %>% dplyr::select(promoter, enhancer) %>% separate(promoter, into = c('chr', 'start', 'end', 'strand'))

tss_df = read.table('./data/raw_data/rnaseq/transcript_loc.bed', header = T, sep = '\t')

overlapped_bed = f_bedtools(f5_promoter[,c('chr', 'start', 'end', 'enhancer')], tss_df[,c('chr', 'start', 'end', 'gene')], fun = 'intersect', paras = '-wao', debug = F)

head(overlapped_bed)

enhancer_gene <- overlapped_bed %>% dplyr::rename( enhancer = V4, gene = V8) %>% dplyr::select(enhancer, gene) %>%
     tidyr::separate(enhancer, into = c('chr', 'start', 'end')) %>% distinct() %>% arrange(gene)

enhancer_gene_filter = enhancer_gene %>% filter(gene != '.')



##Overlap selected features of each gene with the enhancer-gene pairs/ compare with unselected features.
unselected_features = features_stats$merge_control
selected_features = features_stats$features


good_performance_gene <- features_stats$performances %>% filter(performance > 0.05)

selected_bed = f_extract_bed_from_features(selected_features) 
control_bed = f_extract_bed_from_features(unselected_features)

good_gene_list = str_replace(good_performance_gene$gene, '[.][0-9]*', '')
real_overlap = f_overlap_tf_features_with_f5_enhancers(selected_bed, enhancer_gene_filter, good_gene_list)
control_overlap = f_overlap_tf_features_with_f5_enhancers(control_bed, enhancer_gene_filter, good_gene_list)

fisher_table = rbind(table(real_overlap$f5_hit),
table(control_overlap$f5_hit))[,c(2,1)]

colnames(fisher_table) = c('overlap_F5', 'no_overlap')
rownames(fisher_table) = c('selected', 'unselected')


##Output
flog.info('Fisher test results:')
print(fisher_table)
fisher_test=fisher.test(fisher_table)
fisher_test

f_write_paper_results('Fantom 5 Fisher table', fisher_table, output_file)
f_write_paper_results('Fisher p-value', fisher_test$p.value, output_file, scientific = T)
f_write_paper_results('Percentage of selected features overlapped by F5 enhancers', fisher_table[1,1]/sum(fisher_table[1,]), output_file, scientific = T)

######################Fantom 5 part done##############################







#######################Count the number of selected features.####################
if (read_flag == TRUE){
    collect_table = f_feature_enrichment( loc_batch_name, paste0('chr', 1:22), mode_list['TF'], debug = F)
    write.table(collect_table, file = f_p('%s/feature_number_table', results_dir), quote = F, sep = '\t', row.names = F)
}else{
    collect_table = read.table(file = f_p('%s/feature_number_table', results_dir), sep = '\t', header = TRUE)
}
#f_feature_enrichment( loc_batch_name, c('chr22'), mode_list['All'], debug = TRUE)
#f_feature_enrichment( loc_batch_name, c('chr22', 'chr2', 'chr10'), mode_list['TF'])

head(collect_table)
colnames(collect_table) = c('feature', 'selected', 'unselected', 'p-value', 'Odds_ratio')
options(digits=10)
collect_table$location = 'Enhancer'
collect_table$location[grep('promoter', collect_table$feature,)]  = 'Promoter'
collect_table$location[grep('P_E', collect_table$feature,)]  = 'P_E-Interaction'
collect_table$location[grep('TFoverlap', collect_table$feature,)]  = 'TF-interaction'
collect_table$feature = str_replace(collect_table$feature, '(promoter|enhancer|P_E|TFoverlap[.])', '')
table(collect_table$location)


fisher_table <- collect_table %>% group_by(location) %>% dplyr::summarise(selected = sum(selected), unselected = sum(unselected))
f_write_paper_results('Test the depletion of distal features compared with promoter features',fisher_table, output_file)
distal_deplation_test = fisher.test(fisher_table[,2:3])
f_write_paper_results('Fisher test pvalue', f_p('%.4e (0 mean <2.2e-16), coef %s', distal_deplation_test$p.value, distal_deplation_test$estimate ), output_file)

collect_table_write = collect_table
collect_table_write['p-value'] = as.numeric(format(as.numeric(collect_table[,'p-value']), digits = 4, scientific=T))

str(collect_table_write)


collect_table_write['Odds_ratio'] = format(as.numeric(collect_table[,'Odds_ratio']), digits = 2, scientific = T)
head(collect_table_write)

colnames(collect_table_write)

str(collect_table_write)



head(collect_table_write)

collect_table_write = collect_table_write[,c('feature', 'location', "selected",  "unselected",  "p-value", "Odds_ratio")]

collect_table_write$selected = as.numeric(collect_table_write$selected)

collect_table_write_sort = f_sort_by_col(collect_table_write, 'selected', TRUE)
#install.packages('ReporteRs')
#install.packages('gdtools')
#f_table_to_word(collect_table_write)

head(collect_table_write_sort)



plot_data = collect_table_write_sort %>% group_by(feature) %>% dplyr::summarise( tf_selected = sum(selected) ) %>%
    mutate( pos = cumsum(tf_selected) - tf_selected/2 ) %>% arrange(-tf_selected)

head(plot_data, n=20)

plot_data %>% filter(feature == 'PU')

#Pie plot
#ggplot(data=plot_data, aes(x=factor(1), y=tf_selected, fill=factor(feature))) +
#  geom_bar(stat="identity") +
#  geom_text(aes(x= factor(1), y=pos, label = feature), size=10) + guides(fill=FALSE) +  # note y = pos
#  coord_polar(theta = "y")


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





collect_table_tf <- collect_table_write %>% group_by(feature) %>% dplyr::summarise( sum_select = sum(selected) ) %>% arrange(-sum_select)
collect_table_tf$selected = collect_table_tf$sum_select
head(collect_table_tf)


#Features peak count and selected number of features.
collect_table_tf$feature = as.character(collect_table_tf$feature)
test_obj = f_feature_count_correlation(collect_table_tf, peak_count)
f_write_paper_results('Correlation between features peak count and selected number of features.', f_p('P-value: %.4e, coefficient: %s', f_get_exact_pvalue(test_obj), test_obj$estimate), output_file)



##Paper: Top 35 TFs are sufficient
sum(collect_table_tf[1:35,'sum_select'])/sum(collect_table_tf[,'sum_select'])
sum(collect_table_tf[1:5,'sum_select'])/sum(collect_table_tf[,'sum_select'])
f_write_paper_results('First five features account for ', f_p('%s of selected features', sum(collect_table_tf[1:5,'sum_select'])/sum(collect_table_tf[,'sum_select'])), file = output_file)


head(collect_table_tf)

collect_table_tf$feature = factor(collect_table_tf$feature, levels = collect_table_tf$feature )

collect_table_tf_new = collect_table_tf %>% mutate( cum = cumsum(selected) )

head(features_merge)

features_merge$tf_rename = str_replace(features_merge$tf_rename, 'DNASE', 'DHS')

top_feature_data <-  features_merge %>% group_by(tf_rename) %>% dplyr::summarise( TF = length(tf_rename), Gene = length(unique(gene)) ) %>% top_n(n = 10) %>% gather(group, value,-tf_rename)

head(top_feature_data)

tf_ranked <-top_feature_data %>% filter(group == 'TF') %>% arrange(-value) 


top_feature_data$group = factor(top_feature_data$group, levels=c('TF', 'Gene'))

top_feature_data$tf_rename = factor(top_feature_data$tf_rename, tf_ranked$tf_rename)

levels(top_feature_data$group)
levels(top_feature_data$group) <- c('TF count', 'Gene count')

tail(top_feature_data)


length(well_genes)

well_genes

f_write_paper_results('At lease % genes has one DHS', top_feature_data %>% filter(tf_rename == 'DHS', group == 'Gene count' )%>% summarise(value/nrow(well_genes)), output_file)




feature_freq_plot <- 
    ggplot(data=top_feature_data, aes(x=tf_rename, y=value, fill = group)) + #geom_line(aes(y = cum, group =1)) +
  geom_bar(stat="identity", position = 'dodge') + theme_Publication(12) +  theme(axis.text.x = element_text(angle = 90, size = 8,hjust = 1, vjust = 0)) +
    xlab('Features') + ylab('Times selected by TF2Exp models') + theme(legend.position = c(0.8,0.8), legend.direction = 'vertical') +
    scale_fill_discrete("") + theme( axis.text.x = element_text(face="bold", size=10))

feature_freq_plot
ggsave(f_p('%s/feature_freq.tiff', figure_dir), plot = feature_freq_plot, width = 7, height = 4.5)
ggsave(f_p('%s/feature_freq2.tiff', figure_dir), plot = feature_freq_plot, width = 7, height = 7)




##################END of the feature count section########################


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






