
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
        cat('Zero files in A','\n')
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
    batch_name = '54samples_peer'
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

library (plyr)
df <- ldply (str_split(qtl_data$GENE_ID,'[.]'))
qtl_data$gene = df$V1

qtl_bed = qtl_data[,c('CHR_SNP', 'SNPpos','SNPpos','gene')]
qtl_bed$CHR_SNP = paste0('chr', qtl_bed$CHR_SNP)
qtl_bed$snp_id = qtl_data$SNP_ID 


#Read the feature data.
output_file = f_p('%s/output/lasso.prediction.%s.%s.%s.%s.txt', output_dir, chr_str, model_str, ifelse(add_histone,  'addHis', 'rmvHis'), ifelse(add_miRNA, 'addMiR', 'rmvMiR') )
chr_feature_file = f_p('%s.features', output_file)
chr_feature_data  = read.table(chr_feature_file, header = TRUE)
chr_feature_data = chr_feature_data[complete.cases(chr_feature_data),]


head(chr_feature_data)




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
    if (!is.null(overlap_bed)){
        cat(gene_id,'\n')
    }else{
        next
    }

    merge_overlap_bed = rbind(merge_overlap_bed, as.data.frame( overlap_bed))

}


write.table(merge_overlap_bed, file = f_p('%s/output/eqtl_features.txt', output_dir), quote = F)



