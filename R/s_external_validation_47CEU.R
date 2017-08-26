#This script is try to use the data from "Population Variation and Genetic Control of Modular Chromatin Architecture in Humans"
#to validate the selected features in the TF2Exp model.

#Step 1. Get the selected features for each chromosom
source('~/R/s_ggplot2_theme.R')
source('s_summary_fun.R')
source('s_project_funcs.R')
source('s_gene_regression_fun.R')
source('r_feature_analysis_fun.R')
library(optparse)
library(dplyr)

f_associate_two_features<- function(input_data1, input_data2, id1, id2, shared_indivs, debug = F){

    if (debug == TRUE){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        
    }

    stats_data = data.frame( pvalue = rep(1, length(id1)))
    stats_data$coef = 0
    f_ASSERT(length(id1) == length(id2), 'Matched index')
    for (i in 1:length(id1)){
        loc_peak_id = as.character(id1[i])
        loc_gene = as.character(id2[i])
        cor_obj = cor.test( as.numeric(input_data2[loc_gene, shared_indivs]),  as.numeric(input_data1[loc_peak_id, shared_indivs] ), method = 'spearman', exact = FALSE)
        stats_data[i, 'coef'] = cor_obj$estimate
        stats_data[i, 'pvalue'] = cor_obj$p.value
    }
    
    stats_data$fdr = p.adjust(stats_data$pvalue, method = 'BH')
    
    return (stats_data[, c('pvalue', 'coef', 'fdr')])

}



#collection_name = 'peer358corRmdup'
collection_name = 'validation'
loc_batch_name = '358samples_regionkeepLow'



factor_name = 'PU1'
real_factor_name = 'PU.1'


#factor_name = 'CTCF'
#real_factor_name = 'CTCF'


#################Load selected features and model performance#####################
results_dir = f_p('./data/%s/rnaseq/results/%s/', loc_batch_name, collection_name)

##data stored by s_summary_regression_features.R
load(file = f_p('%s/features_stats', results_dir)) 


head(features_stats)

features_merge = subset(features_stats, !is.na(feature_start) & mode == factor_name)
features_merge$tf_rename = str_replace(features_merge$rename, '(enhancer|promoter)', '')
table(features_merge$tf_rename)
features_merge$location = 'promoter'
features_merge$location[ grepl('enhancer', features_merge$rename) ] = 'enhancer'

features_merge$tf_rename[features_merge$tf_rename == 'PU'] = 'PU.1'
features_merge$tf_rename[features_merge$tf_rename == 'USF'] = 'USF.1'

#sample_performance_merge = read.table(file = f_p('%s/sample_performance_merge.%s', results_dir, collection_name), header = T)
sample_performance_merge = read.table(file = f_p('%s/sample_performance_merge', results_dir), header = T)

tf_performance = subset(sample_performance_merge, mode == factor_name)
rownames(tf_performance) = tf_performance$gene
head(tf_performance)

head(features_merge$rename)
sort(table(features_merge$tf_rename))
length(unique(features_merge$gene))


library(tidyr)

tf_feature_data = features_merge %>% filter(tf_rename == real_factor_name) %>%
    group_by(feature_start, feature_end, gene) %>% filter(score == max(score)) %>%
    unite(col = binding_id, chr, feature_start, feature_end, remove = F)

best_features <- features_merge %>% group_by(gene) %>% summarise( best_score = score[which.max(abs(score))],
                                                                 best_feature = tf_rename[which.max(abs(score))], feature_count = length(gene) ) %>% filter(best_feature == real_factor_name)

dim(features_merge)
head(features_merge)
dim(best_features)

#output: tf_feature_data



##################Load deepsea data############################################
deepsea_data = read.table(f_p('/homed/home/shi/expression_var/R/data/358samples_regionkeepLow/output/tf_variation/%s_deepsea.txt', tolower(factor_name)),
                          header = T, sep = '\t', fill = TRUE, row.names = NULL, stringsAsFactors = F)
rownames(deepsea_data) = paste0(deepsea_data$chr, '_', deepsea_data$start, '_', deepsea_data$end)
deepsea_data_sub = deepsea_data[,4:449]
deepsea_data_sub[deepsea_data_sub == '.'] = 0
deepsea_data_sub = data.matrix(deepsea_data_sub)
#Output: deepsea_data_sub









#########################Load PU1 binding data in 47 individuals#######################
binding_data = read.table( f_p('/homed/home/shi/expression_var/R/data/raw_data/CEU47/peer9_%s.txt', factor_name), sep = ' ', fill = TRUE,quote = "", header = TRUE)
rownames(binding_data) = binding_data$PeakID#peak id might have duplicates
#f_ASSERT( length(setdiff(tf_feature_data$binding_id, binding_data$PeakID )) ==0, 'Peak ID match problem!')
#Ouput: binding data






#####Load expression in 47 and 358 individuals##############
tf_individules = grep('NA[0-9].*', colnames(binding_data), value = T)
shared_indivs = tf_individules
gene_expression = read.table('/homed/home/shi/expression_var/R/data/358samples_regionkeepLow/rnaseq/GEUVADIS.Gene.DATA_MATRIX', header = T, sep = ' ', row.names = 1)
expression_47indivs = read.table('/homed/home/shi/expression_var/data/raw_data/CEU47/peer10_gene.txt', header = T, sep = ' ', row.names = 1)

small_indivs = intersect(shared_indivs, colnames(expression_47indivs))
shared_genes=intersect(rownames(gene_expression), rownames(expression_47indivs))

rnaseq_cor = f_associate_two_features(gene_expression, expression_47indivs, id1 = shared_genes, id2 = shared_genes, shared_indivs = small_indivs , debug = F)
rnaseq_cor %>% summarise(sig_count = sum(fdr<0.05, na.rm = T), mean_coef = mean(coef, na.rm = T))
hist(subset(rnaseq_cor, fdr<0.05)$coef)
hist(subset(rnaseq_cor, fdr<1)$coef)

rownames(rnaseq_cor) = shared_genes
stable_genes = shared_genes[rnaseq_cor$fdr<0.05]
shared_indivs = intersect(colnames(gene_expression), tf_individules)





######################Load Peak binding score and reference binding score####################
tf_peak_score = read.table(f_p('/homed/home/shi/expression_var/R/data/raw_data/CEU47/%s_peak_ref.txt', factor_name), header = T )
tf_peak_score$PeakID = with(tf_peak_score, paste(chr, start, end, sep = '_'))
rownames(tf_peak_score) = tf_peak_score$PeakID
f_ASSERT( length(setdiff(tf_feature_data$binding_id, tf_peak_score$PeakID)) == 0, "Peak ID doesn't match")
hist(tf_peak_score$score, breaks = 100, main = 'DeepSEA binding score of PU.1 peaks on ref allele')


peak_variation_stats_44 = rowSums(abs(deepsea_data_sub[tf_peak_score$PeakID, shared_indivs]) > 0.01 )
peak_variation_stats_358 = rowSums(abs(deepsea_data_sub[tf_peak_score$PeakID, ]) > 0.01)
peak_variation_stats_minor_44 = rowSums(abs(deepsea_data_sub[tf_peak_score$PeakID, shared_indivs]) > 1e-4 )

hist(peak_variation_stats_minor_44, breaks =100)

#peak_variation_stats_44 = rowSums(abs(deepsea_data_sub[tf_peak_score$PeakID, shared_indivs]) > 0.02*tf_peak_score$score )

head(deepsea_data_sub)

deepsea_sd = apply(deepsea_data_sub[tf_peak_score$PeakID, shared_indivs], 1, sd)

#hist(deepsea_sd, breaks = 40)


length(peak_variation_stats_44)
hist(peak_variation_stats_44)
#output: peak_variation_stats


##############Deepsea and TF binding##############################


dim(tf_peak_score)

head(binding_data) #normalized 
binding_matrix = read.table(file = f_p('/homed/home/shi/expression_var/R/data/raw_data/CEU47/%s_raw.txt', factor_name), sep = '\t')
rownames(binding_matrix) = binding_matrix$PeakID
binding_matrix$PeakID = NULL

dim(binding_matrix)
head10(binding_matrix)

tf_peak_score$binding_signal = rowMeans(binding_matrix[tf_peak_score$PeakID, ])/(tf_peak_score$end - tf_peak_score$start)
dim(tf_peak_score)
dim(binding_matrix)

head(tf_peak_score)

ggplot(tf_peak_score, aes(binding_signal, score)) + geom_point() + xlab('Mean binding signal') + ylab('Ref DeepSEA score')





#########Overall correlation between deepsea and TF bidning##################
CEU47_deepsea        = deepsea_data_sub[, shared_indivs]
head(CEU47_deepsea)


raw_cor = f_associate_two_features(CEU47_deepsea,  binding_data, binding_data$PeakID, binding_data$PeakID, shared_indivs )
rownames(raw_cor) = binding_data$PeakID


length(unique(binding_data$PeakID))
dim(raw_cor)
dim(binding_data)
raw_cor %>% filter(fdr < 0.05) %>% summarise(neg = sum(coef<0), pos = sum(coef > 0), rate = sum(coef<0)/length(coef) )
raw_cor %>% summarise(fdr_count = sum(fdr < 0.05, na.rm =T), total = length(coef), significant_rate = sum(fdr<0.05, na.rm = T)/length(coef) )



CEU47_deepsea_sub = CEU47_deepsea[which(rowSums(abs(CEU47_deepsea) > 0.01) >= 0.05*ncol(CEU47_deepsea) ),]
CEU47_deepsea_sub = CEU47_deepsea[deepsea_sd[rownames(CEU47_deepsea)] > 0.005,]
dim(CEU47_deepsea_sub)
head(CEU47_deepsea_sub)
dim(CEU47_deepsea)
shared_peaks = intersect(rownames(CEU47_deepsea_sub), binding_data$PeakID)

overall_cor = f_associate_two_features(CEU47_deepsea_sub,  binding_data, shared_peaks, shared_peaks, shared_indivs )
rownames(overall_cor) = shared_peaks

#Overall the wrong correlation between deepsea and TF binding.
overall_cor %>% filter(fdr < 0.05) %>% summarise(neg = sum(coef<0), pos = sum(coef > 0), rate = sum(coef<0)/length(coef) )
overall_cor %>% summarise(fdr_count = sum(fdr < 0.05, na.rm =T), total = length(coef), significant_rate = sum(fdr<0.05, na.rm = T)/length(coef) )


hist(subset(overall_cor, fdr < 0.05)$coef, breaks = 100, main = 'FDR < 0.05')
hist(overall_cor$coef, breaks = 100, main = 'DeepSEA and TF binding correlation in 43k peaks across 44 indiv ')
head(overall_cor)

length(shared_peaks)

length(peak_variation_stats_44)

raw_cor$variant_count = peak_variation_stats_44[rownames(raw_cor)]
raw_cor$variant_count[is.na(raw_cor$variant_count)] =0
raw_cor$deepsea_var = deepsea_sd[rownames(raw_cor)]
raw_cor$ref_score = tf_peak_score[rownames(raw_cor), 'score']

ggplot(raw_cor %>% filter(fdr < 0.05), aes( coef > 0 ,variant_count)) + geom_boxplot() + ylab( 'Count of DeepSEA score diff > 0.01' )
wilcox.test( subset(raw_cor, fdr < 0.05 & coef <0)$variant_count,  subset(raw_cor, fdr < 0.05 & coef >0)$variant_count)

ggplot(raw_cor %>% filter(fdr < 1), aes( coef > 0 ,deepsea_var)) + geom_boxplot()
wilcox.test( subset(raw_cor, fdr < 0.05 & coef <0)$deepsea_var,  subset(overall_cor, fdr < 0.05 & coef >0)$deepsea_var)

head(raw_cor)
raw_cor %>% na.omit() %>% filter(fdr < 0.05)%>% mutate(pos_cor = coef > 0) %>% group_by(pos_cor) %>% dplyr::summarise(mean_sd = mean(deepsea_var))
#ggplot(overall_cor %>% filter(fdr < 0.05), aes( coef > 0 ,ref_score)) + geom_boxplot()

hist(raw_cor$deepsea_var, breaks = 100)

sum(raw_cor$deepsea_var < 0.01, na.rm = T)
sum(raw_cor$deepsea_var < 0.005, na.rm = T)




#############Identify Detectable genes under 40.############
train_samples = setdiff( colnames(gene_expression), shared_indivs)
length(train_samples)

repeat_times = 3
results = data.frame(matrix(0, nrow = nrow(tf_feature_data), ncol = repeat_times))
for ( i in 1:repeat_times){

    selected_samples = sample(train_samples, size = 44, replace = F)
    print(selected_samples)
    test_result = f_associate_two_features( deepsea_data_sub, gene_expression, tf_feature_data$binding_id, tf_feature_data$gene,  shared_indivs, debug = F)
    results[,i] = test_result$fdr
}
detectable_genes = which( rowSums(results < 1) == repeat_times )

library(pwr)
#detectable_coef = pwr.r.test(n = 44, power = 0.8)$r
#detectable_genes = which(abs(test_result$coef) > detectable_coef)
detectable_genes_names = tf_feature_data[detectable_genes, 'gene']$gene
intersect(stable_genes, detectable_genes_names)

dim(tf_performance)

intersect(stable_genes, rownames(subset(tf_performance, performance > 0.05)) )
#output: detectable_gene_names









#############################################
tf_feature_data$stable_genes = tf_feature_data$gene %in% stable_genes
table(tf_feature_data$stable_genes)

tf_feature_data$perf = tf_performance[tf_feature_data$gene,'performance']
tf_feature_data$deepsea_sd = deepsea_sd[tf_feature_data$binding_id]
tf_feature_data$variant_count = peak_variation_stats_44[tf_feature_data$binding_id]


##tf_feature_data_sub = tf_feature_data[detectable_genes,] %>% ungroup() #%>% filter(stable_genes == T)
names(which(peak_variation_stats_44 > 4))

dim(tf_feature_data)
tf_feature_data$alteration_number = peak_variation_stats_minor_44[tf_feature_data$binding_id]
hist(tf_feature_data$alteration_number)





tf_feature_data_sub = tf_feature_data %>% ungroup() %>% filter(binding_id %in% names(which(peak_variation_stats_44 >= 3)))
dim(tf_feature_data_sub)



sd_threshold = quantile(deepsea_sd, 0.90)

tf_feature_data_sub = tf_feature_data %>% ungroup() %>% arrange(-deepsea_sd) %>% filter(deepsea_sd > sd_threshold )#quantile(deepsea_sd, .80))

head(tf_feature_data_sub)
dim(tf_feature_data_sub)
##tf_feature_data_sub = tf_feature_data %>% ungroup() %>% filter(binding_id %in% names(which(deepsea_sd > 0.003)))
#tf_feature_data_sub = tf_feature_data %>% ungroup() %>% filter(alteration_number > 5)

tf_feature_shared_peaks = intersect(rownames(binding_data), intersect(tf_feature_data_sub$binding_id, rownames(deepsea_data_sub)))

tf_feature_data_sub <- tf_feature_data_sub %>% filter(binding_id %in% tf_feature_shared_peaks)
dim(tf_feature_data_sub)

tf_feature_data_sub[,c('tf_deepsea_pvalue', 'tf_deepsea_coef', 'tf_deepsea_fdr')] = f_associate_two_features(binding_data, deepsea_data_sub, tf_feature_data_sub$binding_id, tf_feature_data_sub$binding_id, shared_indivs)

tf_feature_data_sub[,c('gene_deepsea_pvalue', 'gene_deepsea_coef', 'gene_deepsea_fdr')] = f_associate_two_features( deepsea_data_sub, gene_expression, tf_feature_data_sub$binding_id, tf_feature_data_sub$gene,  shared_indivs, debug = F)

tf_feature_data_sub[,c('tf_gene_pvalue', 'tf_gene_coef', 'tf_gene_fdr')] = f_associate_two_features(binding_data, gene_expression, tf_feature_data_sub$binding_id, tf_feature_data_sub$gene, shared_indivs)


tf_feature_data_sub = tf_feature_data_sub %>% ungroup() %>% na.omit()


fdr_threshold = 0.05

tf_feature_data_sub %>% summarise(
                        tf_gene_count = sum(tf_gene_fdr < fdr_threshold, na.rm = T),
                        #tf_gene_pvalue_sig = sum(tf_gene_pvalue < 0.05, na.rm = T),
                        tf_deepsea_count = sum(tf_deepsea_fdr<fdr_threshold, na.rm = T),
                        gene_deepsea_count = sum(gene_deepsea_fdr < fdr_threshold, na.rm = T),
                        total = length(score)
                        )


dim(tf_feature_data_sub)
tf_feature_data_sub %>% ungroup() %>% filter( tf_gene_fdr < fdr_threshold ) %>% mutate(tf_gene_flag = tf_gene_fdr < fdr_threshold,  tf_deepsea_flag = tf_deepsea_fdr < 0.05) %>%
    select(binding_id, gene, stable_genes, score, tf_gene_fdr ,matches('tf.*_coef|_flag|pvalue'), deepsea_sd, variant_count ) %>% as.data.frame()

colnames(tf_feature_data_sub)

dim(tf_feature_data_sub)
tf_feature_data_sub$factor = factor_name
write.table(as.data.frame(tf_feature_data_sub), file = f_p('./data/raw_data/CEU47/%s_validation_results', factor_name), quote = FALSE, row.names =F)


ggplot(tf_feature_data_sub %>% filter( tf_gene_fdr < 1), aes(tf_gene_coef, gene_deepsea_coef, alpha = 1, color = tf_gene_fdr < fdr_threshold )) + geom_point() + geom_abline()

ggplot(tf_feature_data_sub %>% filter( tf_gene_fdr < fdr_threshold), aes(tf_gene_coef, gene_deepsea_coef, alpha = 1, size = tf_deepsea_fdr < fdr_threshold, color = -1 * variant_count )) + geom_point() + geom_abline()

sig_tf_gene_pairs = tf_feature_data_sub %>% filter( tf_gene_fdr < 0.05)
sig_tf_gene_pairs = tf_feature_data_sub



cor.test( sig_tf_gene_pairs$tf_gene_coef, sig_tf_gene_pairs$gene_deepsea_coef, paired = T)

write(x = f_p("%s: %s selected features in %s individuals, %s above variance threshold", factor_name, nrow(tf_feature_data), length(shared_indivs), nrow(tf_feature_data_sub)), file = f_p('./data/raw_data/CEU47/%s_stats.txt', factor_name ))

stop()



#on the TF-Deepsea Results.
f_check_deepsea_tf_correlation <- function(input_data){
    print(input_data %>% ungroup() %>% summarise( tf_deepsea_sig = sum(tf_deepsea_fdr < 0.05), total = length(tf_deepsea_coef), rate = tf_deepsea_sig/total ))
    print(input_data %>% ungroup() %>% filter(tf_deepsea_fdr < 0.05) %>% summarise( tf_deepsea_pos = sum(tf_deepsea_coef > 0), total = length(tf_deepsea_coef), rate = 1 -  tf_deepsea_pos/total ))
    
}

f_check_deepsea_tf_correlation(tf_feature_data_sub)

tf_feature_data_sub$tf_deepsea_pvalue

head(tf_feature_data_sub)

#Similar to the overall correlation

tf_feature_data_sub %>% ungroup() %>% filter( tf_gene_fdr < 0.05) %>% f_check_deepsea_tf_correlation
tf_feature_data_sub %>% ungroup() %>% f_check_deepsea_tf_correlation



#Problem:
#TF_gene and deepsea_gene are different
tf_feature_data_sub %>% filter( tf_gene_fdr < 0.05 | tf_deepsea_fdr< 0.05) %>% dplyr::select(score, tf_gene_coef)

tf_feature_data_sub %>% ungroup() %>% filter( tf_gene_fdr < 0.05, tf_gene_coef > 0, gene_deepsea_coef < 0 ) %>% select( binding_id, matches('_coef') )

ggplot(tf_feature_data_sub %>% filter( tf_gene_fdr < 0.05), aes(tf_gene_coef, gene_deepsea_coef, color = tf_deepsea_fdr < 0.05)) + geom_point() + geom_abline()

ggplot(tf_feature_data_sub %>% filter( tf_gene_fdr < 0.05), aes(gene_deepsea_coef, score, colour = tf_gene_pvalue < 0.05)) + geom_point() + geom_abline()

hist(subset(tf_feature_data_sub, tf_deepsea_fdr < 0.05)$tf_deepsea_coef)


tf_feature_data_sub %>% filter( tf_gene_fdr < 0.05 ) %>% summarise(
                                                             neg_deepsea = sum(tf_deepsea_coef<0),
                                                             deep_cor = sum(tf_deepsea_fdr < 0.05)
                                                         )

stat_data <- tf_feature_data_sub %>% filter(tf_gene_fdr < 0.05) %>% select(tf_deepsea_coef, gene_deepsea_coef, tf_gene_coef, tf_deepsea_fdr)

cor.test(stat_data$tf_gene_coef, stat_data$gene_deepsea_coef)

tf_feature_data_sub$variant_count = peak_variation_stats[tf_feature_data_sub$binding_id]

ggplot(tf_feature_data_sub %>% filter(tf_gene_fdr < 0.05), aes( tf_deepsea_coef < 0 ,variant_count)) + geom_boxplot() + geom_point()

correlated_peaks <- tf_feature_data_sub %>% filter(tf_gene_fdr < 0.05) %>% select(binding_id, gene) 

tf_feature_data_sub$tf_gene_flag = tf_feature_data_sub$tf_gene_fdr < 0.05
tf_feature_data_sub$tf_deepsea_flag = tf_feature_data_sub$tf_deepsea_fdr < 0.05
tf_feature_data_sub$deepsea_coef_flag = tf_feature_data_sub$tf_deepsea_coef > 0

f_make_peak_gene_names <- function(input_data){
    paste0(input_data$binding_id, '_', input_data$gene)
}


annotation_data = data.frame(tf_feature_data_sub)
rownames(annotation_data) = f_make_peak_gene_names(tf_feature_data_sub)

tf_feature_data_sub %>% filter(binding_id == 'chr9_116113247_116113440') %>% as.data.frame

head(tf_feature_data_sub)

tf_feature_data_sub[correlated_peaks$binding_id, c('tf_deepsea_flag', 'tf_deepsea_coef') ]


library(NMF)


aheatmap(deepsea_data_sub[correlated_peaks$binding_id,], annRow = annotation_data[f_make_peak_gene_names(correlated_peaks), c('tf_deepsea_flag', 'deepsea_coef_flag') ],
         filename = '/homed/home/shi/expression_var/R/data/raw_data/CEU47/a.tiff')

tf_feature_data_sub %>% filter(tf_gene_fdr < 0.05,tf_deepsea_coef > 0, tf_deepsea_fdr < 0.05)

target_peak = 'chr1_166938406_166938682'

target_peak = 'chr1_22106751_22107027'
target_peak = 'chr4_25310030_25310200'
qplot( as.numeric(binding_data[target_peak, shared_indivs]), deepsea_data_sub[target_peak,shared_indivs])
tf_peak_score[target_peak,]

class(binding_data)
class(binding_data[1,])
as.vec(binding_data[1,])

sum(deepsea_data_sub[target_peak,shared_indivs] > 0.01)




















