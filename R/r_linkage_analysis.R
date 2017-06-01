source('s_summary_fun.R')
source('r_feature_analysis_fun.R')
library(ggplot2)
source('s_linkage_fun.R')
#ggplot(ld_table, aes(MAF_A, MAF_B)) + geom_point()
source('r_bedtools.R')
source('s_project_funcs.R')





f_get_ld_between_two_feature_set <- function(input_features, selected_genes, chr_str, snp_feature_df, ld_table, debug = FALSE){
    
    if (f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
        browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        
    }
    
    #return_list = list()
    #For each TF feature, idenfity the max impact SNP and number of SNPs.
    tf_features <- input_features %>% filter(gene %in% selected_genes)
    max_snp_table <- tf_features %>% group_by(name) %>% dplyr::summarise(max_snp = snp[which.max(abs(impact_score))], snp_count = length(impact_score)) %>% as.data.frame
    rownames(max_snp_table) = max_snp_table$name
    tf_features$max_snp = tf_features$snp == max_snp_table[tf_features$name, 'max_snp']
    tf_features$snp_count = max_snp_table[tf_features$name, 'snp_count']
    shared_genes = unique( intersect( snp_feature_df$gene, tf_features$gene ) )
    overlapped_snps_in_snp_model = snp_feature_df %>% filter(gene %in% shared_genes)

    #Do the overlapping between SNPs in TF and SNPs in the SNPinTF models
    merge_df_raw= merge(tf_features,  snp_feature_df[,c('gene', 'snp')], by = 'gene', allow.cartesian=TRUE )
    f_ASSERT(length(shared_genes) == length(unique(merge_df_raw$gene)))
    same_hit_regions <- merge_df_raw %>% filter(snp.x == snp.y) %>% mutate(gene_snp_id = paste0(gene, '|SNP.', snp.y))
    head(merge_df_raw)
    head(same_hit_regions)

    
    head(same_hit_regions)
    sel_peaks_with_mulitple_sel_SNPs <- same_hit_regions %>% group_by(name) %>% dplyr::summarise( snp_count = length(snp.y) ) %>% dplyr::summarise( sum(snp_count > 1)/length(snp_count) )
    flog.info( '%.2f  the selected peaks overlap with multiple selected SNPs', sel_peaks_with_mulitple_sel_SNPs[1,1] )

    
    
    snp_with_multiple_tfs=same_hit_regions %>% group_by(snp.y) %>% dplyr::summarise(count = length(snp.y)) %>% dplyr::summarise( sum(count > 1))
    snp_hit_max <- merge_df_raw %>% filter(name %in% same_hit_regions$name, max_snp == TRUE, snp.x == snp.y, snp_count > 1 ) %>% mutate(gene_snp_id = paste0(gene, '|SNP.', snp.y))

    
    return_list = c(snp_with_multiple_tfs[1,1], length(unique(snp_feature_df$snp)),
                    length(unique(same_hit_regions$name)), length(unique(tf_features$name)),
                    sum(max_snp_table$snp_count > 1), nrow(max_snp_table),
                    length(unique(snp_hit_max$snp.y)),
                    length( unique( paste0( same_hit_regions$gene, same_hit_regions$snp.y))),
                    nrow(overlapped_snps_in_snp_model)
                    )

    return_list = unlist(return_list)
    names(return_list) = c('snp_with_multiple_tfs', 'total_snp_num',
                           'tf_overlap_with_snps', 'total_tf_num',
                           'tf_with_multiple_snps', 'total_tf_num_baseds_maxSNP',
                           'max_snp_in_tf', 'snp_selected', 'overlapped_snps_in'

                           
                           )

    f_ASSERT(length(unique(tf_features$name)) == nrow(max_snp_table), 'TF features should be same')
    #setdiff(unique(tf_features$name) , unique(merge_df_raw$name))
    
    cat( as.numeric(snp_with_multiple_tfs[1,1])/length(unique(snp_feature_df$snp)), 'of ', length(unique(snp_feature_df$snp)) , ' key SNPs overlapped with multipe key TF regions \n')
    flog.info('TF model feature regions %.3f of %s are overlapped with selected SNP in SNPinTF models overlapped ,
                  %.3f peaks with multiple SNPs,
                  %.3f of overlapped SNPs total are max with the TF bound peaks',
              length(unique(same_hit_regions$name))/length(unique(tf_features$name)),
              length(unique(tf_features$name)),
              sum(max_snp_table$snp_count > 1)/nrow(max_snp_table), #max_snp_table is based on each TF selected features, snp_count is the number of overlapped SNPs
              length(unique(snp_hit_max$snp.y))/length(unique(snp_feature_df$name)))
    
    flog.info('%.3f of %s selected SNPs overlapped with key TF features', length( unique( paste0( same_hit_regions$gene, same_hit_regions$snp.y)))/nrow(overlapped_snps_in_snp_model),
                  nrow(overlapped_snps_in_snp_model)             )

    #Check the LD score of other snps to the selected SNPs
    merge_df <- merge_df_raw %>% filter(name %in% same_hit_regions$name)

    merge_df$index1 = paste0(merge_df$snp.x, ':', merge_df$snp.y)
    merge_df$index2 = paste0(merge_df$snp.y, ':', merge_df$snp.x)
    index_fit1 = merge_df$index1 %in% row.names(ld_table)
    index_fit2 = merge_df$index2 %in% row.names(ld_table)
    
    merge_df$R2 = 0
    merge_df$R2[index_fit1] = ld_table[merge_df$index1[index_fit1],'R2']
    merge_df$R2[index_fit2] = ld_table[merge_df$index2[index_fit2],'R2']
    cat('Associated with other SNPs:', sum(merge_df$R2 > 0), mean(merge_df$R2[merge_df$R2 > 0]), '\n')

    #merge_df %>% filter( R2 > 0) %>% dplyr::summarize(n = unique(gene))


    #For the selected SNPs overlaped with the TF regions, compare the contribution of closely linked SNP and selected SNPs.
    tf_features$id = paste0(tf_features$name, '|', tf_features$snp )
    rownames(tf_features) = tf_features$id

    same_hit_regions$id = with(same_hit_regions, paste0(name, '|', snp.x))

    head(same_hit_regions)
    
    snp_pairs = merge(tf_features[ !(tf_features$id %in% same_hit_regions$id), ], same_hit_regions[,c( 'snp.y', 'name', 'id')], by = 'name', allow.cartesian=TRUE )
    snp_pairs_b = f_get_ld_R2_based_on_snp_paris(snp_pairs, 'snp', 'snp.y', ld_table) %>% filter(R2 > 0.9)

    plot_data = data.frame(SNP = tf_features[snp_pairs_b$id.y,'impact'], ld_SNP = tf_features[snp_pairs_b$id.x, 'impact']) 

    p <- ggplot(plot_data, aes(SNP, ld_SNP)) + geom_point() + geom_abline() + xlab('Variant impact on the TF')

    print(p)

    dim(snp_pairs_b)
    #return (merge_df)

    return (return_list)
}







#####################Compare features selected by SNP model and TF models.############################################

mode_list = modes_list[['peer358corRmdup']]

SNPinTF =mode_list['SNPinTF']
batch_name = '358samples_regionkeepLow'

collection_name = 'peer358corRmdup'



chr_str_list = paste0('chr', c(22:1))
#chr_str_list = c('chr1', 'chr10', 'chr15', 'chr22')

chr_str = 'chr22'
data_stats = data.frame()

for (chr_str in chr_str_list){

    cat('=============', chr_str,'================= \n')
    
    ld_table = read.table(f_p('./data/raw_data/wgs/1kg/ld_results.%s.ld', chr_str), header = T)
    ld_table$id = paste0(ld_table$SNP_A,  ':', ld_table$SNP_B)
    ld_table = ld_table[!duplicated(ld_table$id),]
    rownames(ld_table) = ld_table$id

    maf_table_raw = read.table(f_p('./data/raw_data/wgs/1kg/%s.freq.frq', chr_str), header = T)

    head(maf_table_raw)
    maf_table = maf_table_raw[ !duplicated(maf_table_raw$SNP), ]
    rownames(maf_table) = maf_table$SNP

    ld_table$MAF_A = maf_table[ld_table$SNP_A, 'MAF' ]
    ld_table$MAF_B = maf_table[ld_table$SNP_B, 'MAF']

    head(ld_table)

    ##return_list1 = f_merge_selected_features(batch_name, chr_str_list, SNPinTF, debug = F)
    return_list1 = f_summary_regression_results(batch_name, chr_str, SNPinTF, rsync_flag = TRUE, return_features = TRUE)

    good_genes = subset(return_list1$performance, performance > 0.05)$gene
    snp_feature_df = subset(return_list1$features, rename == 'SNP')
    head(snp_feature_df)
    ##All the SNP features without performance fitlering.
    snp_feature_df$gene = str_replace(snp_feature_df$name, '[|].*', '')
    snp_feature_df$snp = str_replace(snp_feature_df$name, '.*SNP.', '')
    head(snp_feature_df)


    ##Get the SNVs from the rareVar model
    TF_model = mode_list['TF']
                                        #return_list_tf = f_merge_selected_features(batch_name, chr_str_list, SNPinTF, debug = F)

    return_list_tf = f_summary_regression_results('358samples_regionkeepLow_snpOnly', chr_str, TF_model, rsync_flag = FALSE, return_features = TRUE)
    return_list_tf$control
    head(return_list_tf$features)


    feature_df = return_list_tf$features

    

    var_loc = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.loc', chr_str), sep = '\t', header = T)
    colnames(var_loc) = c('chr', 'pos', 'name', 'ref', 'alt')
    var_maf = read.table(f_p('./data/raw_data/wgs/1kg/maf/%s.maf', chr_str), sep = '\t', header = T)
    var_impact = read.table( f_p('./data/%s/deep_result/all/chrMerge2/evalue/%s.diff.gz', batch_name, chr_str), header = T, sep = ',')

    var_maf$snp_id = paste0(var_maf$pos, ':', var_maf$ref, ':', var_maf$alt)

    var_maf <- var_maf %>% distinct(snp_id, .keep_all = T)

    head(var_maf)

    var_feature_bed_subset = f_snp_and_impact_in_tf_features(return_list_tf$features, chr_str, var_loc, var_impact, var_maf, debug = F)

    var_feature_bed_subset$impact_score = abs( var_feature_bed_subset$impact)


    tf_good_perf <- return_list_tf$performance %>% filter(rsq > 0.05)

    unique(tf_good_perf$gene)


    head(var_feature_bed_subset)

    var_feature_bed_subset_rmdup <-  var_feature_bed_subset %>% distinct(tf_start, tf_end, gene, snp, .keep_all = T )


    shared_predictable_genes = intersect(tf_good_perf$gene, good_genes)


    #For the selected SNPs overlaped with the TF regions, compare the contribution of closely linked SNP and selected SNPs.
    chr_data =    f_get_ld_between_two_feature_set(var_feature_bed_subset_rmdup, shared_predictable_genes, chr_str, snp_feature_df, ld_table, debug =F)

   
    data_stats = rbind(data_stats, c(chr_str, chr_data))
    
    table(var_feature_bed_subset$chr)
    select <- dplyr::select
                                        #Some of the selected SNPs are the top imacted variants in the peaks.
                                        #Some of them are not, but may be linked to the top positions in the peak.
                                        #This function compares the linkage between the selected SNP (not in the top positions) and the rest SNPs.
    f_compare_ld_with_top_impact_positions_for_selected_snp_and_rest_snps(var_feature_bed_subset_rmdup, tf_good_perf$gene, chr_str, snp_feature_df, ld_table, debug =T) 
}


colnames(data_stats) = c('chr', names(chr_data))

summary_results = as.list(colSums(data.matrix(data_stats[,2:ncol(data_stats)])))
names(summary_results) = colnames(data_stats[,2:ncol(data_stats)])

class(summary_results)

summary_results

summary_results$snp_with_multiple_tfs


results_dir = f_p('./data/%s/rnaseq/results/%s/', batch_name, collection_name)
linkage_output_file = f_p('%s/linkage_results.txt', results_dir)
f_write_paper_results('=====linkage analysis======', data = date(), file = linkage_output_file)



cat( as.numeric(summary_results$snp_with_multiple_tfs/summary_results$total_snp_num), 'of ', summary_results$total_snp_num , ' key SNPs overlapped with multipe key TF regions \n',
    file = linkage_output_file, append = T)

names(summary_results)

cat(
    f_p('TF model feature regions %.3f of %s are overlapped with selected SNP in SNPinTF models overlapped ,
                  %.3f peaks with multiple SNPs,
                  %.3f of overlapped SNPs total are max with the TF bound peaks
          ',
    summary_results$tf_overlap_with_snps/summary_results$total_tf_num,
    summary_results$total_tf_num,
    summary_results$tf_with_multiple_snps/summary_results$total_tf_num_baseds_maxSNP,
    summary_results$max_snp_in_tf/summary_results$total_snp_num),
    file = linkage_output_file, append = T)


cat(
f_p('%.3f of %s selected SNPs overlapped with key TF features\n', summary_results$snp_selected/summary_results$total_snp_num,
                  summary_results$total_snp_num), file = linkage_output_file, append = T)


stop()







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











