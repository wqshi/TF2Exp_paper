
#This script combine the results from s_regression_for_one_gene.R

#Sync the data from clustdell.


setwd('~/expression_var/R/')
source('~/R/s_function.R', chdir = TRUE)
source('s_gene_regression_fun.R')
source('s_summary_fun.R')
library(doMC)
library(stringr)
#install.packages('doMC')

library(GenomicRanges)
library("optparse")

option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="dataset file name", metavar="character"),
    make_option(c("--chr"),      type="character", default='chr22', help="chromosome name", metavar="character")    
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$batch_name)){
    #batch_name = '54samples_evalue'
    #batch_name = '462samples_sailfish_quantile'
    #batch_name = '462samples_snyder_norm'
    batch_name = '462samples_quantile'
    #batch_name = '462samples_quantile_rmNA'
    chr_str = 'chr22'
    
}else{
    batch_name = opt$batch_name
    chr_str = opt$chr
}



output_dir = f_p("./data/%s/rnaseq/%s/", batch_name, chr_str)

#Rsync the features files back to the clustdell
rsync_cmd = f_p("rsync -rav  --include '*enet*' --exclude '*' shi@clustdell.cmmt.ubc.ca:/home/shi/expression_var/data/%s/rnaseq/%s/   ~/expression_var/data/%s/rnaseq/%s/", batch_name, chr_str,batch_name, chr_str)
system(rsync_cmd)



features_df=f_merge_data(output_dir, '.*enet.feature')
performance_df = f_merge_data(output_dir, '.*enet$')
head(performance_df)

head(features_df)

head(performance_df)
cat('Mean performance:',mean(performance_df$performance), '\n')
hist(performance_df$performance)

library(dplyr)
good_predictions <- performance_df %>% filter(performance > 0.4)

cat(nrow(good_predictions), 'out of', nrow(performance_df), 'have good predictions(Rsquared > 0.4)','\n')


write.table(performance_df, f_p('%s/performance', output_dir), quote =FALSE, sep = '\t', row.names= FALSE)


stop()
##########################
###Get the bed files for high accuracy and low accuracy genes.
#######################

performance_threshold = 0.1

head(performance_df)

all.entrezgene = read.table('./data/raw_data/rnaseq/all.ensemble.genes',sep='\t', header = TRUE)
row.names(all.entrezgene) = all.entrezgene$ensembl_transcript_id

head(all.entrezgene)
duplicated_rows = duplicated(all.entrezgene$ensembl_gene_id)
all.entrezgene = all.entrezgene[!duplicated_rows,]
row.names(all.entrezgene) = all.entrezgene$ensembl_gene_id

#Change to the promoter regions.
head(all.entrezgene)
rna_seq = all.entrezgene[, c('chromosome_name','transcript_start', 'transcript_end', 'strand')]
rna_seq = rna_seq[complete.cases(rna_seq),]
head(rna_seq)
colnames(rna_seq) = c('chr','start','end', 'strand')
tmp_data = rna_seq[,c('chr','start','end', 'strand')]
forward_strand = rna_seq$strand == 1
reverse_strand  = rna_seq$strand == -1
cat('Forward strand:', sum(forward_strand), 'Reverse ', sum(!forward_strand) ,'\n')

rna_seq[forward_strand,'start'] = tmp_data[forward_strand, 'start'] - 2000
rna_seq[forward_strand,'end'] = tmp_data[forward_strand, 'start'] + 2000
rna_seq[reverse_strand,'start'] = tmp_data[reverse_strand, 'end'] - 2000
rna_seq[reverse_strand,'end'] = tmp_data[reverse_strand, 'end'] + 2000

library(plyr)
gene_ids=ldply(str_split(performance_df$gene, '[.]'))
performance_df[,c('chr','start','end', 'strand')] = rna_seq[gene_ids$V1,]
performance_df = performance_df[complete.cases(performance_df),]
head(performance_df)

performance_threshold = 0.1
table(performance_df$performance > performance_threshold)
performance_df$chr = paste0('chr', performance_df$chr)

high_group = subset(performance_df, performance >= performance_threshold )
low_group = subset(performance_df, performance < performance_threshold )


write.table(high_group[, c('chr','start','end', 'strand')], f_p('%s/high.bed', output_dir), quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
write.table( performance_df[, c('chr','start','end', 'strand')], f_p('%s/low.bed', output_dir), quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)










