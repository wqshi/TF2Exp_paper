setwd('~/expression_var/R/')
source('~/R/s_function.R', chdir = TRUE)
source('s_gene_regression_fun.R')
source('s_summary_fun.R')
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
    batch_name = '462samples_snyder_norm'
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

performance_df %>% filter(performance > 0.4)

write.table(performance_df, f_p('%s/performance', output_dir), quote =FALSE, sep = '\t', row.names= FALSE)
