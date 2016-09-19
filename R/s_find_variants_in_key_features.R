

chr_str = 'chr22'
batch_name = '445samples_sailfish'
source('s_summary_fun.R')
source('r_feature_analysis_fun.R', chdir = T)
maf_table_raw = read.table('./data/raw_data/wgs/1kg/freq_stat.frq', header = T)
library(dplyr)


TF_model = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.TF_other.info.normCor'
return_list_tf = f_summary_regression_results('445samples_sailfish', 'chr22', TF_model, rsync_flag = FALSE, return_features = TRUE)


var_loc = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.loc', chr_str), sep = '\t', header = T)
colnames(var_loc) = c('chr', 'pos', 'name', 'ref', 'alt')
var_maf = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.maf', chr_str), sep = '\t', header = T)
var_impact = read.table( f_p('./data/%s/deep_result/all/chrMerge/evalue/%s.diff.gz', batch_name, chr_str), header = T, sep = ',')




var_feature_bed_subset = f_snp_and_impact_in_tf_features(return_list_tf$features, chr_str, var_loc, var_impact, var_maf)


head(var_feature_bed_subset)


head(var_feature_bed_subset)
var_feature_bed_subset$variant_type = 'common_snp'
var_feature_bed_subset$variant_type[var_feature_bed_subset$maf < 0.05] = 'rare_var' 

head(var_feature_bed_subset)

dim(var_feature_bed_subset)

colnames(var_feature_bed_subset)

var_feature_bed_subset = as.data.frame(var_feature_bed_subset)
var_feature_bed_concise = var_feature_bed_subset[,c('tf_start', 'tf_end', 'gene', 'tf',  'impact', 'maf', 'variant_type', 'snp_id')]

setdiff(c('tf_start', 'tf_end', 'gene', 'tf', 'snp_id', 'impact', 'maf','variant_type'), colnames(var_feature_bed_subset))

head(var_feature_bed_concise, n = 20)

var_feature_bed_concise %>% group_by(variant_type, tf) %>% dplyr::summarize(mean_imp = mean(impact), var_imp = var(impact), snp_count = length(impact), rare_count = length(impact))


variants_stats <- var_feature_bed_concise %>% group_by(gene, tf) %>% dplyr::summarize(max_imp = max(abs(impact)), index = variant_type[which.max(abs(impact))])

flog.info('Table: the variantions in the selected TF regions')
print(table(var_feature_bed_concise$variant_type))

flog.info('Table: max impact variations')
print(table(variants_stats$index))

flog.info('Table: rare vs SNP in the whole population')
print(table(var_maf$MAF > 0.05))




