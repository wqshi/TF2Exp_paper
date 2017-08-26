source('s_gene_regression_fun.R')
library(tidyr)
library(dplyr)
source('s_gene_data_class.R')
source('s_summary_fun.R')
f_my_cor_df_rows <- function(row1, row2){
    cor( unlist( row1), unlist(row2), 'pairwise.complete.obs')
}

#batch_name_A = '445samples_region'
#collection_name = 'tradR2keepZero'
batch_name_A = '358samples_regionkeepLow'
batch_name_B = f_p('%s_snpOnly',batch_name_A)

collection_name = 'peer358'

keepZero_list = modes_list[[collection_name]]

chr_str = 'chr22'


f_compare_correlation_of_top_features <- function(batch_name_A, batch_name_B, return_list_A, return_list_B, chr_str, input_gene, sample_info, ref_binding_score, debug = F){
    test_gene <- GENE(data = data.frame(), gene_name = input_gene, chr_str = chr_str, batch_name = batch_name_A)
    test_gene$read_data()
    test_gene2 <- GENE(data = data.frame(), gene_name = input_gene, chr_str = chr_str, batch_name = batch_name_B)
    test_gene2$read_data()

    if(f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
    
    }
    sample_cols = test_gene$get_samples()

    better_features <- return_list_B$features %>%
        separate(name, into = c('gene', 'feature'), sep = '[|]') %>% mutate( TF = str_replace(feature, '(promoter|enhancer)[.]', '') ) %>%
        filter(gene == input_gene, feature != '(Intercept)') %>% arrange(desc(abs(score)))

    #print(better_features[1:3,])
    feature_data1 <- test_gene$data %>% filter(feature_end == better_features[1,'feature_end'], feature_start == better_features[1,'feature_start'] )
    feature_data2 <- test_gene2$data %>% filter(feature_end == better_features[1,'feature_end'], feature_start == better_features[1,'feature_start'] )

    tf_name = str_replace(better_features[1, 'TF'], '[.][0-9]+', replacement = '')
    deep_cols = grep(tf_name, colnames(ref_binding_score), value = TRUE, ignore.case = T)
    #print(ref_binding_score[,c(deep_cols,'pos')] %>% filter(pos >=better_features[1,'feature_start'], pos <=  better_features[1,'feature_end']))
    
    merge_feature = rbind(feature_data1[1,sample_cols], feature_data2[1,sample_cols], test_gene$data[1, sample_cols])
    merge_feature[is.na(merge_feature)] = 0

    diff_cols=colnames(merge_feature)[merge_feature[1,] != merge_feature[2,]]
    similar_cols=colnames(merge_feature)[merge_feature[1,] == merge_feature[2,]]
    print(table(sample_info[diff_cols,'pop']))
    cat('Diff individuals', length(diff_cols), 'Similar individuals', length(similar_cols), '\n')
    #print(merge_feature[, diff_cols])

    

    
    cat(
        'A batch feature cor overall:',f_my_cor_df_rows(feature_data1[1, sample_cols], merge_feature[3, sample_cols]),
        'B batch feature cor overall', f_my_cor_df_rows(feature_data2[1, sample_cols], merge_feature[3, sample_cols]), '\n',
        'A diff cor',f_my_cor_df_rows(merge_feature[1, diff_cols], merge_feature[3, diff_cols]),
        'B diff cor', f_my_cor_df_rows(merge_feature[2, diff_cols], merge_feature[3, diff_cols]),
        'Similar cor',f_my_cor_df_rows(merge_feature[2, similar_cols], merge_feature[3, similar_cols]), '\n'
    )

    merge_feature_t = t(merge_feature)
    A_all = data.frame(feature = merge_feature_t[sample_cols, 1], exp = merge_feature_t[sample_cols, 3], type = 'A_all')
    B_all = data.frame(feature = merge_feature_t[sample_cols, 2], exp = merge_feature_t[sample_cols, 3], type = 'B_all')
    A_diff = data.frame(feature = merge_feature_t[diff_cols, 1], exp = merge_feature_t[diff_cols, 3], type = 'A_diff')
    B_diff = data.frame(feature = merge_feature_t[diff_cols, 2], exp = merge_feature_t[diff_cols, 3], type = 'B_diff')
    similar = data.frame(feature = merge_feature_t[similar_cols, 2], exp = merge_feature_t[similar_cols,3], type = 'similar')

    combine_data = rbind(A_all, B_all, A_diff, B_diff, similar)
    colnames(combine_data) = c('feature', 'exp', 'type')
    p<-ggplot(combine_data, aes(feature, exp)) + geom_point( alpha = 0.5 ) + geom_smooth() + facet_wrap(~type)
    
}


return_list_A = f_summary_regression_results(batch_name_A, chr_str, keepZero_list['TFsnpMatch'], rsync_flag = FALSE, return_features = TRUE)
return_list_B = f_summary_regression_results(batch_name_B, 'chr22', keepZero_list['TF'], rsync_flag = FALSE, return_features = TRUE)

head(return_list_B$performance)

shared_genes = intersect(return_list_A$performance$gene, return_list_B$performance$gene)

perf_diff = data.frame(gene = shared_genes, perf_A = return_list_A$performance[shared_genes, 'performance'], perf_B= return_list_B$performance[shared_genes, 'performance'])
perf_diff$diff = perf_diff$perf_A - perf_diff$perf_B

perf_diff = f_sort_by_col(perf_diff, 'diff', decr_flag = TRUE)

head(perf_diff)
sample_info = read.table(f_p('./data/462samples/chr_vcf_files/integrated_call_samples_v3.20130502.ALL.panel'), header = TRUE, row.names = 1)

ref_binding_score = read.table('./data/445samples_regionkeepLow/deep_result/all/chrMergeTF/chr22.ref.gz', header = T, sep = ',')
head(ref_binding_score)

head(perf_diff)
i = 3
for( i in 1:10){
    input_gene = perf_diff[i, 'gene']
    print('\n\n')
    cat('=============', input_gene, 'Performance diff', perf_diff[i, 'diff'] ,'=============', '\n')
    try( result <-
    f_compare_correlation_of_top_features(batch_name_A, batch_name_B, return_list_A, return_list_B, chr_str, input_gene, sample_info, ref_binding_score, debug = F)
    )
}





head(perf_diff)
minus_genes = perf_diff %>% arrange(diff)
head(minus_genes)
rownames(minus_genes) = minus_genes$gene

for(input_gene in minus_genes[1:10, 'gene']){
    print('\n\n')
    cat('=============', input_gene, 'Performance diff', minus_genes[input_gene, 'diff'] ,'=============', '\n')
    try( result <-
    f_compare_correlation_of_top_features(batch_name_A, batch_name_B, return_list_A, return_list_B, chr_str, input_gene, sample_info, ref_binding_score, debug = F)
    )
}


stop()

#21212487    21212867 ENSG00000099940.7 promoter.CHD2 -0.2653038
target_tf = 'CHD2'
target_col=grep(f_p('GM12878.*%s', target_tf), colnames(ref_data), value = T)[1]

ref_data=read.table(f_p('./debug/%s/infile.vcf.out.ref', target_tf), header = T, sep = ',')
alt_data=read.table(f_p('./debug/%s/infile.vcf.out.alt', target_tf), header = T, sep = ',')
diff_data=read.table(f_p('./debug/%s/infile.vcf.out.diff', target_tf), header = T, sep = ',')
ref_data[, c('chr', 'pos', target_col)] %>% filter(pos > 21212487, pos < 21212867) %>% arrange(pos)
alt_data[, c('chr', 'pos', target_col)] %>% filter(pos > 21212487, pos < 21212867) %>% arrange(pos)
diff_data[, c('chr', 'pos', target_col)] %>% filter(pos > 21212487, pos < 21212867) %>% arrange(pos)


perf_diff = f_sort_by_col(perf_diff, 'diff', decr_flag = FALSE)
target_tf = 'YY1'
target_col=grep(f_p('GM12878.*%s', target_tf), colnames(ref_data), value = T)[1]

ref_data=read.table(f_p('./debug/%s/infile.vcf.out.ref', target_tf), header = T, sep = ',')
alt_data=read.table(f_p('./debug/%s/infile.vcf.out.alt', target_tf), header = T, sep = ',')
diff_data=read.table(f_p('./debug/%s/infile.vcf.out.diff', target_tf), header = T, sep = ',')

bed_data=read.table(f_p('./debug/%s/yy1.bed.out', target_tf), header = T, sep = ',')
low_score_data=read.table(f_p('./debug/%s/region_yy1.bed', target_tf), header = T, sep = ',') 

bg_col=grep(f_p('GM12878.*%s', 'ATF2'), colnames(ref_data), value = T)[1]

head10(ref_data)

ref_data[, c('chr', 'pos', target_col)]
alt_data[, c('chr', 'pos', target_col)]
diff_data[, c('chr', 'pos', target_col)]

bed_data[, c('chr','start','end',target_col)]
low_score_data[, c('chr','start',target_col)]

hist(as.numeric(bed_data[,target_col]))
hist(as.numeric(bed_data[,bg_col]))
