
source('~/R/s_function.R', chdir = TRUE)
library("optparse")

option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="dataset file name", metavar="character"),
    make_option(c("--remove_histone"), type="character", default=NULL, help="output file name [default= %default]", metavar="character"),
    make_option(c("--test"), type="character", default=NULL, help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$batch_name)){
    batch_name = '462samples_sailfish_quantile'
    remove_histone = TRUE
    test_flag = TRUE
}else{
    batch_name = opt$batch_name
    remove_histone = opt$remove_histone
    test_flag = opt$test == 'TRUE'
}

chr_str = 'chr22'
model_str = 'all'


output_dir = f_p('./data/%s/', batch_name)
figure_dir = f_p('%s/output/figures/',output_dir)

#Read the expression data for each gene
expression_file = paste0(output_dir, 'rnaseq/gene_regulatory_region_feature_profile.', chr_str )
expression_data = read.table(expression_file, header = TRUE, na.strings = '.')
sample_cols = grep('(NA|HG)[0-9]+', colnames(expression_data),value = T)
expression_data = expression_data[,c('chr', 'start', 'end', 'gene','feature','feature_start', 'feature_end', 'hic_fragment_id' ,'type', sample_cols )]

head(expression_data)
tail(expression_data[,1:10])
dim(expression_data)
cat('Number of transcripts', length(unique(expression_data$gene)), '\n')
cat('Number regulatory regions and empty fragment ids', '\n')
table(expression_data$type, expression_data$hic_fragment_id == '.')
table(expression_data$type, expression_data$feature)
table(unlist(table(expression_data$gene)))
colnames(expression_data)





#Get the expressed transcripts

#write.table(sample_cols, './data/output/sample.list', sep = '\t', quote = F, col.names = F)

rna_data = subset(expression_data, feature == 'RNASEQ')
rownames(rna_data) = rna_data$gene
expressed_transcripts = rownames(rna_data)[rowSums(rna_data[, sample_cols]) > 0]
cat('The number of expressed transcripts:', length(expressed_transcripts), '\n')
head(rna_data)
#str(rna_data)



#Get the genes with at least one elements
gene_with_regulatory_elements = as.data.frame.matrix(table(expression_data$gene, expression_data$feature))
head(gene_with_regulatory_elements)
str(gene_with_regulatory_elements, max.level = 1)
regulated_genes_names = rownames(gene_with_regulatory_elements)[ rowSums(gene_with_regulatory_elements) > 20]
cat('The number of regulated transcripts:', length(regulated_genes_names), '\n')

#Get the inverstigeed genes
genes_names = intersect(expressed_transcripts, regulated_genes_names)
cat('The number of investigated transcripts:', length(genes_names), '\n')

print(genes_names[1:20])

prediction_performance = data.frame(gene = 'mean', performance = '0', stringsAsFactors = F)



sample_info = read.table('/homed/home/shi/expression_var/data/raw_data/wgs/1kg/chr_vcf_files/integrated_call_samples_v3.20130502.ALL.panel', header = TRUE, row.names = 1)
head(sample_info)

dim(sample_info[sample_cols,])

library(dplyr)
dim(expression_data)

colnames(expression_data[,1:15])
head(expression_data[,1:10])

expression_data[is.na(expression_data)] = 0
sample_cols = grep('(NA|HG)[0-9]+', colnames(expression_data), value = TRUE)

library(doMC)
registerDoMC(cores = 7)
zero_cols = (rowSums(expression_data[, sample_cols]==0) >50)
sum(zero_cols)
head(expression_data[zero_cols, sample_cols])
expression_data_subset = expression_data[!zero_cols,]

head(expression_data_subset[,1:10])

gene_stats <- expression_data_subset %>% group_by(gene) %>% summarise(nrow = length(chr), pol2 = sum(feature == 'POL2') , promoter = sum(type == 'promoter'), features = length(unique(feature))  ,hic_id = length(unique(hic_fragment_id)))
gene_stats = as.data.frame(gene_stats)
row.names(gene_stats) = gene_stats$gene
gene_stats$gene = NULL



#Get the regression accuracy according to the pol2 overlapping.
library(caret)



zero_var_genes = nearZeroVar(t(rna_data[, sample_cols]), uniqueCut = 50)
filtered_genes = rownames(rna_data)[zero_var_genes]

remove_histone = FALSE
output_file = f_p('%s/output/lasso.prediction.%s.%s.%s.addMiR.txt', output_dir, chr_str, model_str, ifelse(remove_histone, 'rmvHis', 'addHis') )
performance_data = read.table(file = output_file, sep = '\t')
row.names(performance_data) = performance_data$gene
cat("Mean R-square:",  mean(performance_data[complete.cases(performance_data),]$performance), '\n')
dim(performance_data)
length(filtered_genes)
performance_data[filtered_genes,]



gene_stats[, 'performance'] = as.numeric(performance_data[ rownames(gene_stats) ,'performance'])
gene_stats$overlap_pol2 = as.numeric(gene_stats$pol2) > 0

tail(gene_stats)

table(gene_stats$pol2 == 0)


gene_stats[complete.cases(gene_stats),] %>% group_by(overlap_pol2) %>% summarise(mean =length(performance))

mean_expression = apply(rna_data[,sample_cols],1,mean)
gene_stats$gene_expression = mean_expression[rownames(gene_stats)]
gene_stats$gene_expression_var = apply(rna_data[,sample_cols],1,var)[rownames(gene_stats)]


###Check the correlation between performance and expression mean and var.
colnames(gene_stats)
ggplot(gene_stats, aes(x = performance, y = gene_expression )) + geom_point() + ylim(c(0,1000))
ggsave(f_p('%s/performance_expression.png', figure_dir))


ggplot(gene_stats, aes(x = performance, y = gene_expression_var )) + geom_point() + ylim(c(0,1000))
ggsave(f_p('%s/performance_expression_var.png', figure_dir))

ggplot(gene_stats, aes(x = performance, y = pol2 )) + geom_point()
ggsave(f_p('%s/performance_pol2_num.png', figure_dir))

head(gene_stats)


###




for(loc_gene in rownames(rna_data)){
    qplot(x = as.numeric(rna_data[loc_gene, sample_cols])) + geom_density() + xlab(loc_gene)
    ggsave(f_p('%s/expression/%s.png', figure_dir, loc_gene))
}











