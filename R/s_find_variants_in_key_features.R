source('s_project_funcs.R')
chr_str = 'chr22'
#batch_name = '445samples_region'
source('s_summary_fun.R')
source('r_feature_analysis_fun.R', chdir = T)
maf_table_raw = read.table('./data/raw_data/wgs/1kg/freq_stat.frq', header = T)
library(dplyr)


#TF_model = 'rm.histone_model.cv.glmnet_rm.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.TF_other.info.tradR2keepZeronoInteract'
batch_name = '358samples_regionkeepLow'
TF_model = "rm.histone_model.cv.glmnet_rm.penalty_rm.YRI_population.None_new.batch.358samples.snyder.norm_batch.mode.TF_other.info.tradR2keepZero"


return_list_tf = f_summary_regression_results(batch_name, 'chr22', TF_model, rsync_flag = TRUE, return_features = TRUE)
sum(return_list_tf$performance$performance > 0.05, na.rm = T)


#Read loc of variants, MAF
var_loc = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.loc', chr_str), sep = '\t', header = T)
colnames(var_loc) = c('chr', 'pos', 'name', 'ref', 'alt')
var_maf = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.maf', chr_str), sep = '\t', header = T)
var_impact = read.table( f_p('./data/%s/deep_result/all/chrMerge2/evalue/%s.diff.gz', batch_name, chr_str), header = T, sep = ',')
var_ref = read.table( f_p('./data/%s/deep_result/all/chrMerge2/evalue/%s.ref.gz', batch_name, chr_str), header = T, sep = ',')


#Assign the SNPs impact and ref score to the feature region.
var_feature_bed_subset = f_snp_and_impact_in_tf_features(return_list_tf$features, chr_str, var_loc, var_impact, var_maf, debug = F)
dim(var_feature_bed_subset)
head(var_feature_bed_subset)
var_feature_bed_ref = f_snp_and_impact_in_tf_features(return_list_tf$features, chr_str, var_loc, var_ref, var_maf, debug = F)
var_feature_bed_subset$ref = var_feature_bed_ref$impact
f_ASSERT(all(var_feature_bed_ref$snp == var_feature_bed_subset$snp), 'Mis-match')

var_feature_bed_subset$alt_ratio = with(var_feature_bed_subset, abs(impact)/ref)
var_feature_bed_subset$sign_ratio = with(var_feature_bed_subset, impact/ref)

var_feature_bed_subset$variant_type = 'common_snp'
var_feature_bed_subset$variant_type[var_feature_bed_subset$maf < 0.05] = 'rare_var' 



hist(var_feature_bed_subset$ref, main = 'Reference binding score', xlab = 'ref')
hist(var_feature_bed_subset$impact)

var_feature_bed_subset = as.data.frame(var_feature_bed_subset)

head(var_feature_bed_subset)


##Get the average overlapped variants in the region, For the rare variant part in the manuscript.
table(var_feature_bed_subset$variant_type)/length(unique(var_feature_bed_subset$name))

var_feature_bed_subset_sort = f_sort_by_col(var_feature_bed_subset, 'impact')
head(var_feature_bed_subset_sort)
tail(var_feature_bed_subset_sort)

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




##################Confirm the variation impact calcualtion in one individual for one gene##########
#In this test, I use the genotype of each individual and the variants overlaped in a TF region to recalculate the TF alteration score.

target_gene = 'ENSG00000025708.8'
library(tidyr)
source('s_gene_data_class.R')

##Confirm the variation impact in one TF region is correct.

good_predictions = return_list_tf$performance %>% filter(performance > 0.05)

for( target_gene in good_predictions$gene[1:10] ){    
    f_test_preprocess_for_one_gene(target_gene, chr_str, batch_name, return_list_tf, var_feature_bed_subset, debug = F)
}
###################Test done#####################################################







##Test the HiC fragments are complete.
##Confirm one TF peaks overlapped with HiC fragments (in distal and poximal hic fragments) are complete for each gene.

##Get all the promoter interaction hic-ids from the interaction data.
full_hic_contact = read.table('./data/raw_data/hic/high_interaction_hic.addChr.chr22.txt')
colnames(full_hic_contact) = c('chr', 'score', 'hic_fragment_id', 'pair')
head(full_hic_contact)

full_hic_loc = read.table('./data/raw_data/hic/fragmentStartEnd.addChr.addMissing.txt')
colnames(full_hic_loc) = c('chr', 'start', 'end', 'name')
head(full_hic_loc)

for( target_gene in good_predictions$gene[1:10]){
    f_test_one_tf_complete_in_gene_processed_data(target_gene, batch_name, full_hic_loc, full_hic_contact, return_list_tf, debug = T)
}



if (FALSE){
#####Find the missing part###########
missing_hic_fragment = c('17577374')

missing_hic_fragment %in% features_in_processed$pair
missing_hic_fragment %in% features_in_processed$hic_fragment_id
full_hic_contact %>% filter(pair %in% missing_hic_fragment, score >= 0.4) %>% filter(hic_fragment_id %in% promoter_hic_ids)
full_hic_contact %>% filter(hic_fragment_id %in% missing_hic_fragment, score >= 0.4) %>% filter(pair %in% promoter_hic_ids)

#HET impact for each individual this is not necessary
indiv_het_impact = read.table(f_p('%s/data/%s/deep_result/all/%s/evalue/%s_het/%s.diff', project_dir, batch_name, chr_str, batch_name, individual_id ),
                              sep = ',', header = T)

indiv_het_impact$snp_id = paste(indiv_het_impact$pos, indiv_het_impact$ref, indiv_het_impact$alt, sep = ':')
rownames(indiv_het_impact) = indiv_het_impact$snp_id
deepsea_cols = grep( f_p('.*%s', key_tfs), colnames(indiv_het_impact), value = T, ignore.case = T)

head10(indiv_het_impact)
deepsea_cols
sum(indiv_het_impact[individual_snps_in_tf, deepsea_cols ])

}




