
#install.packages('bnlearn')
source('~/R/s_function.R', chdir = TRUE)
library(doMC)

library(stringr)

library("optparse")

option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="dataset file name", metavar="character"),
    make_option(c("--add_histone"), type="character", default='TRUE', help="output file name [default= %default]", metavar="character"),
    make_option(c("--add_miRNA"), type="character", default='TRUE', help="Add miRNA or not", metavar="character"),
    make_option(c("--test"), type="character", default=NULL, help="output file name [default= %default]", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

f_overlap_two_bed_dataframe <- function(bed_A, bed_B){
    if(nrow(bed_A) ==0) {
        cat('Zero files in A','\n')
        return (NULL)
     }
     if(nrow(bed_B) ==0) {
        cat('Zero files in B','\n')
        return (NULL)
     }
    file_a = tempfile()
    file_b = tempfile()
    write.table(bed_B, file = file_b, sep='\t', col.names = F, quote = F, row.names = F)
    write.table(bed_A, file = file_a, sep='\t', col.names = F, quote = F, row.names = F)
    
    library(data.table)
    cmd <- f_p("bedtools intersect -wao -a %s -b %s", file_a, file_b)
    df <- fread(cmd)

    unlink(file_a)
    unlink(file_b)

    return (df)

}



if (is.null(opt$batch_name)){
    #batch_name = '54samples_peer'
    batch_name = '462samples_quantile_rmNA'
    add_histone = TRUE
    add_miRNA = TRUE
    test_flag = TRUE
}else{
    batch_name = opt$batch_name
    add_histone = opt$add_histone == 'TRUE'
    add_miRNA = opt$add_miRNA == 'TRUE'
    test_flag = opt$test == 'TRUE'
    cat('Test flag :', test_flag, '\n')
}

chr_str = 'chr22'
model_str = 'all'




#Read the qtl data
output_dir = f_p('./data/%s/', batch_name)
qtl_file        = './data/raw_data/qtl/GEUV/YRI89.gene.cis.FDR5.all.rs137.txt'
qtl_data = read.table(qtl_file, header = TRUE)
head(qtl_data)
qtl_data$SNPpos = as.integer(qtl_data$SNPpos)

library (plyr)
df <- ldply (str_split(qtl_data$GENE_ID,'[.]'))
qtl_data$gene = df$V1

qtl_bed = qtl_data[,c('CHR_SNP', 'SNPpos','SNPpos','gene')]
qtl_data$SNPpos
qtl_bed$SNPpos = as.character(qtl_bed$SNPpos)
qtl_bed$SNPpos.1 = as.character(qtl_bed$SNPpos.1 + 1)

qtl_bed$CHR_SNP = paste0('chr', qtl_bed$CHR_SNP)
qtl_bed$snp_id = qtl_data$SNP_ID 


#Read the feature data.
mode_name = 'add.histone_rm.miRNA_rm.TF.exp_model.glmnet_rm.permutation_rm.TF.exp.only_add.predict.TF'
chr_feature_file = f_p('%s/rnaseq/%s/%s/features', output_dir, chr_str, mode_name)
chr_feature_data  = read.table(chr_feature_file, header = TRUE)
chr_feature_data = chr_feature_data[complete.cases(chr_feature_data),]

head(chr_feature_data)

performance_file = f_p('%s/rnaseq/%s/%s/performance', output_dir, chr_str, mode_name)
performance = read.table(performance_file, header = TRUE)
head(performance)

gene_name_full_list =ldply( str_split(as.character(chr_feature_data$name),'[|]') )$V1
gene_list =  unique(ldply(str_split(gene_name_full_list, '[.]'))$V1)

head(gene_list)


gene_id = gene_list[1]
gene_id = 'ENSG00000138964'


#Merge the 

merge_overlap_bed = data.frame()

for (gene_id in gene_list){
    cat('=============', gene_id, '===========')
    bed_A = subset( qtl_bed, grepl(gene_id, qtl_bed$gene ))
    bed_B = subset(chr_feature_data, grepl(gene_id, chr_feature_data$name))
    overlap_bed =f_overlap_two_bed_dataframe(bed_A, bed_B)
    dim(bed_A)
    if (!is.null(overlap_bed)){
        cat(gene_id,'\n')
    }else{
        next
    }
    merge_overlap_bed = rbind(merge_overlap_bed, as.data.frame( overlap_bed))

}

str(merge_overlap_bed)

colnames(merge_overlap_bed) = c('chr', 'start', 'end', 'gene', 'snp_id','chr_B','start_B', 'end_B', 'name_B', 'score_B', 'overlap_length')

cat(length(unique(merge_overlap_bed$gene)), 'Out of ', length(gene_list), 'genes', '\n')
cat(sum(merge_overlap_bed$overlap_length), 'Out of ', nrow(merge_overlap_bed), 'rsnps', '\n')
library(dplyr)
eqtl_stats <- merge_overlap_bed %>% group_by(gene) %>% dplyr::summarise(overlap = sum(overlap_length), eqtl_num = length(overlap_length))
eqtl_stats
rownames(performance) = str_replace(performance$gene, '[.][0-9]*', '')

eqtl_stats$performance = performance[equl_stats$gene, 'performance']

cat('Mean performance for the gene with eQTL:',mean(eqtl_stats$performance), '\n')
cat('Mean performance for the gene with eQTL:',mean(eqtl_stats[ eqtl_stats$overlap > 0, ]$performance), '\n')

write.table(merge_overlap_bed, file = f_p('%s/output/eqtl_features.txt', output_dir), quote = F)
















