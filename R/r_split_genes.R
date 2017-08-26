setwd('~/expression_var/R/')
source('~/R/s_function.R', chdir = TRUE)
source('s_gene_regression_fun.R')
source('s_project_funcs.R')
library(stringr)
library("optparse")
library(dplyr)
option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="dataset file name", metavar="character"),
    make_option(c("--chr"),      type="character", default='chr22', help="chromosome", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$batch_name)){
    #batch_name = '54samples_evalue'
    #batch_name = '462samples_sailfish_quantile'
    #batch_name = '445samples_diff'
    #batch_name = '462samples_quantile_rmNA'    
    #chr_str = 'chr22'

    #batch_name = '358samples_regionkeepLow'
    #chr_str = 'chr10'

    batch_name = 'test'
    chr_str = 'chr22'
    
}else{
    batch_name = opt$batch_name
    chr_str = opt$chr
}



output_dir = f_p('./data/%s/', batch_name)
expression_file = paste0(output_dir, 'rnaseq/gene_regulatory_region_feature_profile.', chr_str )
SNP_hic_id = read.table( f_p('./data/raw_data/wgs/1kg/additive_dir/assign_hic_id_to_SNP.%s', chr_str), header = TRUE, na.strings = '')
#head(SNP_hic_id)
gene_enhancer_table = read.table( f_p( '%s/rnaseq/gene_regulatory_fragment.%s', output_dir, chr_str), header = TRUE, na.strings = '')
expression_data = read.table(expression_file, header = TRUE, na.strings = '.')
#head(expression_data)
#tail(expression_data[,1:20])

#count the number of fragments each gene has.
expression_data %>% group_by(gene) %>% filter(type == 'promoter') %>%
dplyr::summarise(frag_count = length(unique(hic_fragment_id))) %>% filter(frag_count > 3)

cat('Number regulatory regions and empty fragment ids', '\n')
table(expression_data$type, expression_data$hic_fragment_id == '.')


#table(expression_data$type, expression_data$feature)

colnames(expression_data)


#Get the expressed transcripts
sample_cols = grep('(NA|HG)[0-9]+', colnames(expression_data),value = T)

cat('The number of samples:', length(sample_cols), '\n')
expression_data = expression_data[,c('chr', 'start', 'end', 'gene','feature','feature_start', 'feature_end', 'hic_fragment_id' ,'type', 'cor', 'pair' ,sample_cols )]
#write.table(sample_cols, './data/output/sample.list', sep = '\t', quote = F, col.names = F)
rna_data = subset(expression_data, feature == 'RNASEQ')
rownames(rna_data) = rna_data$gene

#expressed_transcripts = rownames(rna_data)[rowSums(rna_data[, sample_cols]) > 0] #This is wrong when I use the stadardized normalization.
expressed_transcripts = rownames(rna_data)

cat('The number of expressed transcripts:', length(expressed_transcripts), '\n')

#head(rna_data)
#str(rna_data)
#Get the genes with at least one elements
gene_with_regulatory_elements = as.data.frame.matrix(table(expression_data$gene, expression_data$feature))
dim(gene_with_regulatory_elements)
print(gene_with_regulatory_elements)
regulated_genes_names = rownames(gene_with_regulatory_elements)#[ rowSums(gene_with_regulatory_elements) > 10]



cat('The number of regulated transcripts:', length(regulated_genes_names), '\n')

#Get the inverstigeed genes
genes_names = intersect(expressed_transcripts, regulated_genes_names)
cat('The number of investigated transcripts:', length(genes_names), '\n')

length(genes_names)

non_sample_cols = setdiff(colnames(expression_data), sample_cols)
print(non_sample_cols)

#Read the TF expression data.
tf_gene_id = read.table('./data/raw_data/rnaseq/tf_ensemble_id.txt', header = T, stringsAsFactors = FALSE)
rownames(tf_gene_id) = tf_gene_id$tf
expression_data$feature_tf =''
str(expression_data)
head(expression_data[,1:10])

expression_data$feature_tf = (tf_gene_id[as.character(expression_data$feature),'external_name'])


gene_data_dir=f_p( '%s/rnaseq/%s/', output_dir, chr_str)
dir.create(gene_data_dir)
all.entrezgene = f_get_all.entrezgene('./')
gtX = f_get_genotype_matrix(chr_str)


for (i in 1:length(genes_names)){
    
    transcript_id = genes_names[i]
    cat('\n','============',i, transcript_id, '==============','\n')
    transcript_data = expression_data[expression_data$gene == transcript_id,]
    gene_snps=f_get_genotype_matrix_for_gene(transcript_id, debug = F)
    merge_data = f_merge_two_dfs(transcript_data, gene_snps)
    table(merge_data$type)
    write.table(merge_data, file = f_p('%s/%s.txt', gene_data_dir, transcript_id), sep = '\t', quote = FALSE, row.names =FALSE)
}

head(expression_data)













