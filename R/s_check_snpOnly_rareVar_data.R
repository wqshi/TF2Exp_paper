

source('s_project_funcs.R')
source('s_gene_data_class.R')
library(futile.logger)
library(plyr)
library(dplyr)
source('s_gene_regression_fun.R')


f_gene_stats <- function(input_gene){
    cat('nrows', nrow(input_gene$data), '\n')
    

    rmdup_data <- input_gene$data %>% dplyr::arrange(gene, -cor ) %>% distinct(feature, feature_start, feature_end, gene, .keep_all = T) %>% as.data.frame

    print(rmdup_data %>% filter(type != 'SNP') %>% group_by(feature) %>% dplyr::summarise(count = length(feature)))


    promoter_data <- rmdup_data %>% filter(type == 'promoter')
    
    print(rmdup_data %>% filter(feature == 'CHD2', type == 'promoter') %>% select(feature, feature_start, feature_end, cor))
    print(table(rmdup_data$type))

    print(promoter_data %>% group_by(feature) %>% dplyr::summarise(count = length(feature)))
}

target_batch = '358samples_regionkeepLow'
gene_name = 'ENSG00000249222.1'
chr_str= 'chr22'
gene_name = "ENSG00000100376.7"


data_dir =  f_p('./data/%s/rnaseq/%s/', target_batch, chr_str)

gene_files = list.files(data_dir, pattern = '.txt')

selected_files = sample(gene_files, size = 50)


for (gene_file in selected_files){

    
    gene_name = str_replace(gene_file, '.txt', '')
    flog.info('==%s==', gene_name)
    
    test_gene <- GENE(data = data.frame(), gene_name = gene_name, chr_str = chr_str, batch_name = target_batch)
    test_gene$read_data()

    rare_gene <- GENE(data = data.frame(), gene_name = gene_name, chr_str = chr_str, batch_name = f_p( '%s_rareVar', target_batch))
    rare_gene$read_data(debug = F)

    snp_gene <- GENE(data = data.frame(), gene_name = gene_name, chr_str = chr_str, batch_name = f_p( '%s_snpOnly', target_batch))
    snp_gene$read_data(debug = F)

                                        #f_gene_stats(test_gene)
                                        #f_gene_stats(rare_gene)
    test_gene$data$feature_name = with(test_gene$data, paste(feature, feature_start, feature_end, sep = '_'))
    rare_gene$data$feature_name = with(rare_gene$data, paste(feature, feature_start, feature_end, sep = '_'))
    snp_gene$data$feature_name = with(snp_gene$data, paste(feature, feature_start, feature_end, sep = '_'))


    f_ASSERT( length(setdiff(rare_gene$data$feature_name, test_gene$data$feature_name)) == 0, 'Rare < Normal')
    f_ASSERT( length(setdiff(snp_gene$data$feature_name, test_gene$data$feature_name)) == 0, 'SNP < Normal')

    f_ASSERT( length(setdiff(test_gene$data$feature_name, snp_gene$data$feature_name)) > 0, 'Normal > SNP')

}






