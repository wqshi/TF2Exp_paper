
#This script combine the results from s_regression_for_one_gene.R

#Sync the data from clustdell.


setwd('~/expression_var/R/')
source('~/R/s_function.R', chdir = TRUE)
source('s_gene_regression_fun.R')
source('s_summary_fun.R')
library(stringr)
#install.packages('doMC')
library(dplyr)
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
    #batch_name = '462samples_quantile'
    batch_name = '462samples_quantile_rmNA'
    #chr_str = 'chr22'
    chr_str = 'chrTF'
    
}else{
    batch_name = opt$batch_name
    chr_str = opt$chr
}

mode_name = 'glmnet'

f_summary_regression_results <- function(batch_name, chr_str, mode_name){
    output_dir = f_p("./data/%s/rnaseq/%s/%s/", batch_name, chr_str, mode_name)

    #Rsync the features files back to the clustdell
    rsync_cmd = f_p("rsync -rav  --include '*enet*' --exclude '*' shi@clustdell.cmmt.ubc.ca:/home/shi/expression_var/data/%s/rnaseq/%s/%s/   ~/expression_var/data/%s/rnaseq/%s/%s/", batch_name, chr_str, mode_name, batch_name, chr_str, mode_name)
    system(rsync_cmd)

    
    
    features_df=f_merge_data(output_dir, '.*enet.feature', quite = FALSE)
    
    performance_df = f_merge_data(output_dir, '.*enet$')
    head(performance_df)

    head(features_df)

    head(performance_df)
    cat('Mean performance:',mean(performance_df$performance), '\n')
    hist(performance_df$performance)


    good_predictions <- performance_df %>% filter(performance > 0.5)
    
    cat(nrow(good_predictions), 'out of', nrow(performance_df), 'have good predictions(Rsquared > 0.5)','\n')


    write.table(performance_df, f_p('%s/performance', output_dir), quote =FALSE, sep = '\t', row.names= FALSE)
    return (performance_df)
}

batch_list = c('462samples_quantile', '462samples_quantile_rmNA', '462samples_var_quantile', '462samples_snyder_norm')
for (batch_name in batch_list){
    cat('======', batch_name, '========\n')
    #performance = f_summary_regression_results(batch_name, chr_str, mode_name)
    Sys.sleep(1)
}


mode_list = c('add.histone_add.miRNA_add.TF.exp_model.glmnet_add.permutation',
    'add.histone_add.miRNA_add.TF.exp_model.glmnet_rm.permutation',
    'add.histone_add.miRNA_rm.TF.exp_model.glmnet_rm.permutation',
      'add.histone_rm.miRNA_add.TF.exp_model.glmnet_rm.permutation',
    'add.histone_add.miRNA_add.TF.exp_model.glmnet_add.permutation_add.TF.exp.only'
    )

mode_list =c('add.histone_add.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only',
    'add.histone_add.miRNA_rm.TF.exp_model.glmnet_rm.permutation_add.TF.exp.only')



for (loc_mode in mode_list){
    cat('======', loc_mode, '========\n')
    performance = f_summary_regression_results('462samples_quantile_rmNA', chr_str, loc_mode)
    Sys.sleep(1)
}



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










