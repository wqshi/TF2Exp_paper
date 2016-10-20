source('s_project_funcs.R')
chr_str = 'chr22'
batch_name = '445samples_region'
source('s_summary_fun.R')
source('r_feature_analysis_fun.R', chdir = T)
maf_table_raw = read.table('./data/raw_data/wgs/1kg/freq_stat.frq', header = T)
library(dplyr)


TF_model = 'rm.histone_model.cv.glmnet_rm.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.TF_other.info.tradR2keepZeronoInteract'
return_list_tf = f_summary_regression_results(batch_name, 'chr22', TF_model, rsync_flag = TRUE, return_features = TRUE)

sum(return_list_tf$performance$performance > 0.05, na.rm = T)


#Read loc of variants, MAF
var_loc = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.loc', chr_str), sep = '\t', header = T)
colnames(var_loc) = c('chr', 'pos', 'name', 'ref', 'alt')
var_maf = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.maf', chr_str), sep = '\t', header = T)
var_impact = read.table( f_p('./data/%s/deep_result/all/chrMerge2/evalue/%s.diff.gz', batch_name, chr_str), header = T, sep = ',')
var_ref = read.table( f_p('./data/%s/deep_result/all/chrMerge2/evalue/%s.ref.gz', batch_name, chr_str), header = T, sep = ',')

head10(var_impact)
var_impact %>% filter(pos == 46645148)
var_loc %>% filter(pos == 46645148)


var_feature_bed_subset = f_snp_and_impact_in_tf_features(return_list_tf$features, chr_str, var_loc, var_impact, var_maf, debug = F)
dim(var_feature_bed_subset)
head(var_feature_bed_subset)
var_feature_bed_ref = f_snp_and_impact_in_tf_features(return_list_tf$features, chr_str, var_loc, var_ref, var_maf, debug = F)
var_feature_bed_subset$ref = var_feature_bed_ref$impact

f_ASSERT(all(var_feature_bed_ref$snp == var_feature_bed_subset$snp), 'Mis-match')
var_feature_bed_subset$alt_ratio = with(var_feature_bed_subset, abs(impact)/ref)
var_feature_bed_subset$sign_ratio = with(var_feature_bed_subset, impact/ref)
as.data.frame(var_feature_bed_subset)[which(var_feature_bed_subset$alt_ratio > 0.5),c('ref', 'impact', 'alt_ratio','sign_ratio')]

as.data.frame(var_feature_bed_subset)[var_feature_bed_subset$ref < 0.1,c('ref', 'impact', 'alt_ratio','sign_ratio')]

head(var_feature_bed_subset, n = 20)

hist(var_feature_bed_subset$ref, breaks = 100)

head(var_feature_bed_subset)

head(var_feature_bed_ref)

hist(var_feature_bed_ref$impact)
hist(var_feature_bed_subset$impact)

head(var_feature_bed_subset)
var_feature_bed_subset$variant_type = 'common_snp'
var_feature_bed_subset$variant_type[var_feature_bed_subset$maf < 0.05] = 'rare_var' 

head(var_feature_bed_subset)

dim(var_feature_bed_subset)

colnames(var_feature_bed_subset)

var_feature_bed_subset = as.data.frame(var_feature_bed_subset)


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



######################################
###Confirm the variation impact calcualtion in one individual for one gene##########
target_gene = 'ENSG00000025708.8'
library(tidyr)
source('s_gene_data_class.R')

##Confirm the variation impact in one TF region is correct.

good_predictions = return_list_tf$performance %>% filter(performance > 0.05)

for( target_gene in good_predictions$gene[1:10] ){    
    f_test_preprocess_for_one_gene(target_gene, chr_str, batch_name, return_list_tf, var_feature_bed_subset, debug = FALSE)
}


##Test the HiC fragments are complete.
##Confirm one TF peaks overlapped with HiC fragments are complete in the processed data.
##Check the HiC pairs.
##Get all the promoter interaction hic-ids from the interaction data.
full_hic_contact = read.table('./data/raw_data/hic/high_interaction_hic.addChr.chr22.txt')
colnames(full_hic_contact) = c('chr', 'score', 'hic_fragment_id', 'pair')
head(full_hic_contact)

full_hic_loc = read.table('./data/raw_data/hic/fragmentStartEnd.addChr.addMissing.txt')
colnames(full_hic_loc) = c('chr', 'start', 'end', 'name')
head(full_hic_loc)

for( target_gene in good_predictions$gene[1:10]){
    f_test_one_tf_complete_in_gene_processed_data(target_gene, batch_name, full_hic_loc, full_hic_contact, return_list_tf, debug = F)
}


f_test_one_tf_complete_in_gene_processed_data <- function(target_gene, batch_name, full_hic_loc, full_hic_contact, return_list_tf, debug = F){

    if (debug){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
        
    }
    
    ##Read the gene data.
    test_gene <- GENE(data = data.frame(), gene_name = target_gene, chr_str = chr_str, batch_name = batch_name)
    test_gene$read_data()
    test_gene$subset_features_to_snp_contain_region('TF', debug = F)
    test_data_flag=test_gene$read_test_data(test_batch = '800samples', debug = FALSE)
    dim(test_gene$data)
    gene_pos_list=test_gene$hic_within_1Mb(return_pos= TRUE)
    transcript_data = test_gene$data
    rownames(transcript_data) = make.names(paste0(transcript_data$type, '-',transcript_data$feature), unique = T)


    ##Identify the first key features.
    key_features_df <- return_list_tf$features %>% separate(name, into =c('gene', 'feature'), sep ='[|]') %>% filter(gene == target_gene)
    key_features = grep('Intercept',key_features_df$feature,value = T, invert = T)
    if (length(key_features) != 1){
        key_features = key_features[1]
    }
    key_tfs = transcript_data[key_features,'feature']
    key_tfs


    ##Get all the promoter HiC ids of the gene.
    promoter_data = transcript_data %>% filter(type == 'promoter')
    enhancer_data = transcript_data %>% filter(type == 'enhancer') 
    promoter_hic_ids = unique(c(promoter_data$hic_fragment_id, enhancer_data$pair))
    promoter_hic_ids


    ##For the key features, use TF peak to overlap with all the promoter and associated regions
    full_interaction_pair <- full_hic_contact %>% filter(pair %in% promoter_hic_ids, score > 0.4) %>% select(pair, hic_fragment_id)
    full_interaction_hic <- full_hic_contact %>% filter(hic_fragment_id %in% promoter_hic_ids, score > 0.4) %>% select(pair, hic_fragment_id)
    dim(full_interaction_pair)
    dim(full_interaction_hic)
    head(full_interaction_pair)
    head(full_hic_loc)
    full_unique_hic_ids = unique(c(full_interaction_pair$hic_fragment_id, full_interaction_hic$pair))

    full_hic_ids_1M = full_unique_hic_ids[ full_unique_hic_ids >= gene_pos_list[1] & full_unique_hic_ids <= gene_pos_list[2] ]

    interact_hic_bed  <- full_hic_loc %>% filter(name %in% full_hic_ids_1M, chr == 'chr22' ) 


    ##tf_peak_file = list.files('./data/raw_data/tf/encode_peaks/processed/', pattern = f_p('.*%s.*narrowPeak', tolower(key_tfs)) , full.names = T)
    tf_var_file = list.files('./data/445samples_region/output/tf_variation/all/chr22/', pattern = f_p('%s.*txt', tolower(key_tfs)) , full.names = T)
    tf_bed_raw = read.table(tf_var_file, header = T)
    sample_cols = grep('(NA|HG)', colnames(tf_bed_raw), value = T)
    informative_rows = rowSums(tf_bed_raw[, sample_cols] != '.') >= 0.05 * 446

    tf_bed = tf_bed_raw[informative_rows,c(3,1,2,4)]

    overlapped_bed = f_bedtools(tf_bed, interact_hic_bed, fun = 'intersect', paras = '-wb')
    colnames(overlapped_bed) = c('chr', 'feature_start', 'feature_end', 'name', 'chr2', 'hic_start', 'hic_end', 'hic_fragment_id')

    unique_tf_hic=unique(as.data.frame(overlapped_bed)[,c('feature_start', 'feature_end', 'hic_fragment_id')])
    unique_tf_hic$name = with(unique_tf_hic, paste(feature_start, feature_end, hic_fragment_id, sep='_'))

    features_in_processed <- transcript_data %>% filter(feature == key_tfs) %>% select(feature, feature_start, feature_end ,type, cor ,hic_fragment_id, pair)
    unique_feature_regions = unique(features_in_processed[,c('feature_start', 'feature_end' ,'hic_fragment_id')])
    unique_feature_regions$name = with(unique_feature_regions, paste(feature_start, feature_end, hic_fragment_id, sep='_'))

    print('+++++++++++++++++++')
    ##This is possible, because not all the peaks will overlap with an variant.
    if (length(setequal(unique_feature_regions$name, unique_tf_hic$name))){
        
        flog.info('Test past %s', key_tfs)
    }else{
        flog.error('Test failed %s', key_tfs)
        print(setdiff(unique_feature_regions$name, unique_tf_hic$name))
        print(setdiff(unique_tf_hic$name, unique_feature_regions$name))
    }

}



if (FALSE){
#####Find the missing part###########
missing_hic_fragment = c('17577374')

missing_hic_fragment %in% features_in_processed$pair
missing_hic_fragment %in% features_in_processed$hic_fragment_id
full_hic_contact %>% filter(pair %in% missing_hic_fragment, score >= 0.4) %>% filter(hic_fragment_id %in% promoter_hic_ids)
full_hic_contact %>% filter(hic_fragment_id %in% missing_hic_fragment, score >= 0.4) %>% filter(pair %in% promoter_hic_ids)
}



#HET impact for each individual this is not necessary
indiv_het_impact = read.table(f_p('%s/data/%s/deep_result/all/%s/evalue/%s_het/%s.diff', project_dir, batch_name, chr_str, batch_name, individual_id ),
                              sep = ',', header = T)

indiv_het_impact$snp_id = paste(indiv_het_impact$pos, indiv_het_impact$ref, indiv_het_impact$alt, sep = ':')
rownames(indiv_het_impact) = indiv_het_impact$snp_id
deepsea_cols = grep( f_p('.*%s', key_tfs), colnames(indiv_het_impact), value = T, ignore.case = T)

head10(indiv_het_impact)
deepsea_cols
sum(indiv_het_impact[individual_snps_in_tf, deepsea_cols ])
