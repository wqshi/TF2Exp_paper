

chr_str = 'chr22'
source('s_summary_fun.R')
maf_table_raw = read.table('./data/raw_data/wgs/1kg/freq_stat.frq', header = T)
library(dplyr)


TF_model = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.TF_other.info.normCor'
return_list_tf = f_summary_regression_results('445samples_sailfish', 'chr22', TF_model, rsync_flag = FALSE, return_features = TRUE)


var_loc = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.loc', chr_str), sep = '\t', header = T)
colnames(var_loc) = c('chr', 'pos', 'name', 'ref', 'alt')
var_maf = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.maf', chr_str), sep = '\t', header = T)
var_impact = read.table( f_p('./data/%s/deep_result/all/chrMerge/evalue/%s.diff.gz', batch_name, chr_str), header = T, sep = ',')

dim(var_maf)
length(unique(var_maf$name))


head(var_maf)

head(var_loc)
dim(var_maf)


head(var_loc)

var_bed = data.frame(chr = paste0(var_loc$chr), start = var_loc$pos -1, end = var_loc$pos, name = paste0(var_loc$ref, ':', var_loc$alt))

head(var_bed)

f_get_vars_in_key_modle_features <- function(input_features, chr_str, var_bed){
    tf_bed = f_extract_bed_from_features(input_features, chr_str)
    overlapped_bed = f_bedtools(tf_bed, var_bed, fun = 'intersect', paras = '-wao')
    head(overlapped_bed)
    colnames(overlapped_bed) = c('chr', 'tf_start', 'tf_end', 'name', 'chr.2', 'snp_start', 'snp_end', 'snp', 'overlap')
    overlapped_bed$gene = str_replace(overlapped_bed$name, '[|].*', '')
    return (overlapped_bed)
}


batch_name = '445samples_sailfish'
source('r_bedtools.R')
var_feature_bed =f_get_vars_in_key_modle_features(return_list_tf$features, chr_str, var_bed)
var_feature_bed$tf = str_replace(var_feature_bed$name, '.*[|](enhancer|promoter).','')
var_feature_bed$tf = str_replace(var_feature_bed$tf, '[.][0-9]+$', '')
var_feature_bed$tf[var_feature_bed$tf == 'PU'] = 'PU.1'
var_feature_bed$tf[var_feature_bed$tf == 'USF'] = 'USF.1'
var_feature_bed$tf[var_feature_bed$tf == 'EGR'] = 'EGR.1'



rownames(var_impact) = paste0(var_impact$pos, ':', var_impact$ref,':', var_impact$alt)
var_feature_bed$snp_id = paste0(var_feature_bed$snp_end, ':', var_feature_bed$snp)

dim(var_feature_bed)
shared_ids = intersect(var_feature_bed$snp_id, rownames(var_impact))
length(shared_ids)


var_feature_bed_subset = var_feature_bed[var_feature_bed$snp_id %in% shared_ids,]
dim(var_feature_bed_subset)

concise_deepsea_cols = unique(str_replace(colnames(var_impact), 'None.[0-9]+$', 'None'))

length(concise_deepsea_cols)
unique(var_feature_bed_subset$tf)

dim(var_feature_bed_subset)

var_feature_bed_subset$deepsea=ldply(lapply(paste0(var_feature_bed_subset$tf,'.[0-9]{0,1}(None|TNFa)$'), grep, concise_deepsea_cols, value = T, ignore.case = T))


var_feature_bed_subset$impact =0
for (i in 1:nrow(var_feature_bed_subset)){
    var_feature_bed_subset$impact[i] = var_impact[var_feature_bed_subset$snp_id[i] ,var_feature_bed_subset$deepsea[i]]
}


rownames(var_maf) = paste0(var_maf$pos, ':', var_maf$ref, ':', var_maf$alt)

var_feature_bed_subset$maf = var_maf[var_feature_bed_subset$snp_id,'MAF']


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




