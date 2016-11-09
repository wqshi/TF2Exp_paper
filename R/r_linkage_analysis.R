source('s_summary_fun.R')
source('r_feature_analysis_fun.R')
ld_table = read.table('./data/raw_data/wgs/1kg/ld_results.ld', header = T)
row.names(ld_table) = paste0(ld_table$SNP_A,  ':', ld_table$SNP_B)

maf_table_raw = read.table('./data/raw_data/wgs/1kg/freq_stat.frq', header = T)

head(maf_table_raw)
maf_table = maf_table_raw[ !duplicated(maf_table_raw$SNP), ]
rownames(maf_table) = maf_table$SNP


head(ld_table)

ld_table$MAF_A = maf_table[ld_table$SNP_A, 'MAF' ]
ld_table$MAF_B = maf_table[ld_table$SNP_B, 'MAF']


head(ld_table)

library(ggplot2)

#ggplot(ld_table, aes(MAF_A, MAF_B)) + geom_point()

###Get the SNPs from the SNPinTF model:
SNPinTF = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.SNPinTF_other.info.normCor'
return_list1 = f_summary_regression_results('445samples_sailfish', 'chr22', SNPinTF, rsync_flag = FALSE, return_features = TRUE)
head(return_list1$features)
good_genes = subset(return_list1$performance, performance > 0.05)$gene
snp_feature_df = subset(return_list1$features, rename == 'SNP')
head(snp_feature_df)
##All the SNP features without performance fitlering.
snp_feature_df$gene = str_replace(snp_feature_df$name, '[|].*', '')
snp_feature_df$snp = str_replace(snp_feature_df$name, '.*SNP.', '')
head(snp_feature_df)


###Get the SNVs from the rareVar model
TF_model = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.TF_other.info.normCor'
return_list_tf = f_summary_regression_results('445samples_snpOnly', 'chr22', TF_model, rsync_flag = FALSE, return_features = TRUE)
return_list_tf$control
head(return_list_tf$features)


feature_df = return_list_tf$features
chr_str = 'chr22'
source('r_bedtools.R')


batch_name = '445samples_region'

var_loc = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.loc', chr_str), sep = '\t', header = T)
colnames(var_loc) = c('chr', 'pos', 'name', 'ref', 'alt')
var_maf = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.maf', chr_str), sep = '\t', header = T)
var_impact = read.table( f_p('./data/%s/deep_result/all/chrMerge2/evalue/%s.diff.gz', batch_name, chr_str), header = T, sep = ',')

var_feature_bed_subset = f_snp_and_impact_in_tf_features(return_list_tf$features, chr_str, var_loc, var_impact, var_maf)

var_feature_bed_subset$impact_score = abs(var_feature_bed_subset$maf * var_feature_bed_subset$impact)

head(var_feature_bed_subset)





f_get_ld_between_two_feature_set <- function(input_features, selected_genes, chr_str, snp_feature_df, ld_table, debug = FALSE){

    if (f_judge_debug(debug)){
    ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
    }
    tf_features <- input_features %>% filter(gene %in% selected_genes)
    max_snp_table <- tf_features %>% group_by(name) %>% dplyr::summarise(max_snp = snp[which.max(impact_score)]) %>% as.data.frame
    rownames(max_snp_table) = max_snp_table$name
    tf_features$max_snp = tf_features$snp == max_snp_table[tf_features$name, 'max_snp']

    shared_genes = unique( intersect( snp_feature_df$gene, tf_features$gene ) )
    overlapped_snps_in_snp_model = snp_feature_df %>% filter(gene %in% shared_genes)
    
    merge_df_raw= merge(tf_features,  snp_feature_df[,c('gene', 'snp')], by = 'gene', allow.cartesian=TRUE )

    f_ASSERT(length(shared_genes) == length(unique(merge_df_raw$gene)))
    
    same_hit_regions <- merge_df_raw %>% filter(snp.x == snp.y)
    merge_df <- merge_df_raw %>% filter(name %in% same_hit_regions$name)
 
    snp_hit_max <- merge_df %>% filter(max_snp == TRUE, snp.x == snp.y)

    flog.info('Selected SNP in SNPinTF models overlapped TF model feature regions %.3f out of %s, %.3f are max',
              length(unique(same_hit_regions$name))/length(unique(merge_df_raw$name)),length(unique(merge_df_raw$name)),
              length(unique(snp_hit_max$name))/length(unique(merge_df_raw$name)))
    flog.info('%.3f variants out of %s variants', length( unique( paste0( same_hit_regions$gene, same_hit_regions$snp.y)))/nrow(overlapped_snps_in_snp_model),
                  dim(overlapped_snps_in_snp_model)             )

    #Check the LD score of other snps to the selected SNPs
    merge_df$index1 = paste0(merge_df$snp.x, ':', merge_df$snp.y)
    merge_df$index2 = paste0(merge_df$snp.y, ':', merge_df$snp.x)
    index_fit1 = merge_df$index1 %in% row.names(ld_table)
    index_fit2 = merge_df$index2 %in% row.names(ld_table)
    
    merge_df$R2 = 0
    merge_df$R2[index_fit1] = ld_table[merge_df$index1[index_fit1],'R2']
    merge_df$R2[index_fit2] = ld_table[merge_df$index2[index_fit2],'R2']
    cat('Associated with other SNPs:', sum(merge_df$R2 > 0), mean(merge_df$R2[merge_df$R2 > 0]), '\n')

    #merge_df %>% filter( R2 > 0) %>% dplyr::summarize(n = unique(gene))

    #return (merge_df)
}





tf_good_perf <- return_list_tf$performance %>% filter(rsq > 0.05)

unique(tf_good_perf$gene)


f_get_ld_between_two_feature_set(var_feature_bed_subset, tf_good_perf$gene, chr_str, snp_feature_df, ld_table)

fg_linkage = f_get_ld_between_two_feature_set(return_list_tf$features, tf_good_perf$gene, chr_str, snp_feature_df, ld_table)
bg_linkage = f_get_ld_between_two_feature_set(return_list_tf$control, tf_good_perf$gene, chr_str, snp_feature_df, ld_table)
###Calculate the linkage
fg_linkage$type = 'TF'
bg_linkage$type = 'Control'

dim(fg_linkage)
dim(bg_linkage)

plot_data = rbind(fg_linkage, bg_linkage) %>% filter(R2 > 0)

plot_data$linkage = plot_data$R2

ggplot(plot_data, aes(linkage, color = type)) + geom_density() + xlab('Linkage Score') + ggtitle('Variants in selected TF regions \n linked with selected SNPs in SNPinTF model')


plot_data %>% group_by(type) %>% summarise(mean_linkage = mean(linkage), count = length(linkage), TF_snp = length(unique(snp.x)),  linked_snp = length(unique(snp.y)))

fg_linkage %>% summarise(gene = length(unique(gene)), snp = length(unique(snp.y)))

snp_feature_df %>% summarise(gene_count = length(unique(gene)), snp_count = length(snp), unique_snp = length(unique(snp)))


#So most of the selected regions are not overlapped with the SNPinTF key positions.

###Check example####
head(snp_feature_df, n = 8)
tf_good_perf$gene
snp_feature_df %>% group_by(gene) %>% summarise(snp = length(snp))

target_gene = tf_good_perf$gene[3]

tf_features = return_list_tf$features

tf_features[grep(target_gene, tf_features$name),]

snp_feature_df[grep(target_gene, snp_feature_df$gene),]




