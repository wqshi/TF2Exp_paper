
source('~/R/s_function.R', chdir = TRUE)
sailfish_dir = f_p('%s/expression_var/sailfish/', Sys.getenv("HOME"))
batch_name='qsub_47samples'
raw_file_list = list.files(path = f_p('%s/data/%s/gene_quant/', sailfish_dir, batch_name), full.names = TRUE)
expression_files = grep('.*gene_level.sf', raw_file_list, value = TRUE)

i = 1

merge_data = NULL

gene_file = expression_files[1]

library('stringr')

for (gene_file in expression_files){
    cat(gene_file, '\n')
    individual_id =str_replace( basename(gene_file), 'quant.gene_level.sf', '')
    gene_expression_data = read.table(gene_file)
    individual_data = data.frame( TPM = gene_expression_data$V3, gene = gene_expression_data$V1 )
    colnames(individual_data) = c(individual_id, 'gene')
    head(individual_data)
    if (is.null(merge_data)){
        merge_data = individual_data
    }else{
        merge_data=merge(merge_data, individual_data, all = TRUE, by = 'gene')
    }

    cat('After merging, data size:', dim(merge_data),'\n')
}

cat('After merging, data size:', dim(merge_data),'\n')
write.table(merge_data,  f_p('%s/data/%s/gene_quant/GEUVADIS.Gene.DATA_MATRIX', sailfish_dir, batch_name), row.names = FALSE, quote = FALSE)



stop()


####Compare two processes####
exp_370 =read.table(f_p('%s/data/qsub/gene_quant/GEUVADIS.Gene.DATA_MATRIX', sailfish_dir), header = T)
#exp_445 =read.table(f_p('%s/data/qsub_445samples/gene_quant/GEUVADIS.Gene.DATA_MATRIX', sailfish_dir))
exp_445 = merge_data
dim(exp_445)
dim(exp_370)
head(exp_370)
head(exp_445)
cor_list = c()
samples = colnames(exp_370)[-1] # Remove the gene
loc_sample = samples[2]

for (loc_sample in samples){

    cor_list = c(cor_list, cor(x=exp_370[,loc_sample], y = exp_445[, loc_sample], use = 'pairwise.complete.obs'))
    
}

mean(cor_list) #Mean correlation between two approaches 0.9886. Should be fine.

###############Check the single file in the fastq#####
fastq_list=read.table('/homed/home/shi/expression_var/sailfish/data/fastq/fastq.20160718.list')$V1
print(fastq_list)
head(fastq_list)
library(stringr)
file_id=str_replace(fastq_list, '_.*', '')
table(file_id)

#Only the file ERR188406




