################
#Add the chr coordinate using biomart to the GEUVADIS.Transcript.DATA_MATRIX file


#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#install.packages('pbkrtest')
#install.packages('quantreg')
#packageurl <- "http://cran.r-project.org/src/contrib/Archive/car/car_2.0-25.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
#install.packages('caret')
#install.packages('car')
source('~/expression_var/R/s_project_funcs.R')
library(biomaRt)
library(stringr)
library(caret)
source('~/R/s_function.R', chdir = TRUE)

library("optparse")

option_list = list(
    make_option(c("--batch_name"),      type="character", default=NULL, help="dataset file name", metavar="character"),
    make_option(c("--sep"), type="character", default='tab', help="Seperator of the data file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$batch_name)){
    batch_name = '54samples_genebody'
    #batch_name = '462samples_quantile_rmNA'
    #batch_name = '462samples_log_quantile'
    seperator = ' '
}else{
    batch_name = opt$batch_name
    seperator = ifelse( opt$sep == 'tab', '\t', ' ') 
}

#batch_name = '54samples_gene'
#seperator = ' '

#batch_name = '462samples'
#seperator = '\t'
cat(batch_name, '\n')
batch_output_dir = f_p('./data/%s/', batch_name)


read_flag = FALSE
transcript_flag = FALSE
if (read_flag){
    # define biomart object
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice"  )
    mart_attributes=(listAttributes(mart))
    mart_attributes[grep('end_position', mart_attributes$name, ignore.case = T),]
    # query biomart
    #all.entrezgene <- unique( getBM(attributes = c("ensembl_gene_id", 'ensembl_transcript_id', 'chromosome_name','transcript_start', 'transcript_end', 'strand', 'entrezgene'), values = "*", mart = mart) )
    #all.entrezgene <- unique( getBM(attributes = c("ensembl_gene_id", 'ensembl_transcript_id', 'chromosome_name','transcript_start', 'transcript_end', 'strand', 'uniprot_genename', 'external_gene_name', 'wikigene_name'), values = "*", mart = mart) )
    all.entrezgene <- unique( getBM(attributes = c("ensembl_gene_id", 'chromosome_name','start_position', 'end_position', 'strand', 'uniprot_genename', 'external_gene_name'), values = "*", mart = mart) )
    all.entrezgene$ensembl_transcript_id = all.entrezgene$ensembl_gene_id
    colnames(all.entrezgene) = c("ensembl_gene_id", 'chromosome_name','transcript_start', 'transcript_end', 'strand', 'uniprot_genename', 'external_gene_name', 'ensembl_transcript_id')
    all.entrezgene = all.entrezgene[!duplicated(all.entrezgene$ensembl_gene_id),]
    table(all.entrezgene$chromosome_name)
    all.entrezgene = all.entrezgene[grep('MT|_|GL',all.entrezgene$chromosome_name, invert = T),]

    row.names(all.entrezgene) = all.entrezgene$ensembl_transcript_id    
    dim(all.entrezgene)
    head(all.entrezgene,n=20)
    all.entrezgene[all.entrezgene==''] = '.'
    head(all.entrezgene, n=20)
    
    write.table(all.entrezgene, file ='./data/raw_data/rnaseq/all.ensemble.genes.gene_start',sep='\t', quote = FALSE, row.names = FALSE )
}else{
    cat('Skip query biomart data','\n')
    all.entrezgene = read.table('./data/raw_data/rnaseq/all.ensemble.genes.gene_start',sep='\t', header = TRUE, quote = "")
    row.names(all.entrezgene) = all.entrezgene$ensembl_transcript_id
}

#attributes = listAttributes(mart)
#head(attributes)
#grep('entrez', attributes$name)
#attributes[55,]


#chr_22 = subset(all.entrezgene, chromosome_name == '22')
#head(chr_22)


#write.table(unique(chr_22$ensembl_gene_id), file = 'gene.id.txt', quote = F, sep = '\t', row.names = FALSE, col.names = FALSE)

#sum(is.na(chr_22$entrezgene))
#dim()

head(all.entrezgene)
table(all.entrezgene$strand)
#table(all.entrezgene$chromosome_name)

colnames(all.entrezgene)
"ENST00000000233" %in% all.entrezgene$ensembl_transcript_id
head(all.entrezgene)

rna_seq_raw = read.csv(file = f_p('%s/rnaseq/GEUVADIS.Gene.DATA_MATRIX', batch_output_dir) , sep=seperator, header =TRUE)
cat('The size of the raw rna_seq data is:', dim(rna_seq_raw), '\n')

#colnames(rna_seq_raw)

sample_cols = grep('(NA|HG)[0-9]+', colnames(rna_seq_raw),value = T)
print(sample_cols[1:10])
rna_expression_values = as.matrix(rna_seq_raw[, sample_cols])
dim(rna_expression_values)
#head(rna_expression_values)
near_zero = caret::nearZeroVar(t(rna_expression_values))

if(length(near_zero) == 0){
    rna_seq = rna_seq_raw
}else{
    rna_seq = rna_seq_raw[-(near_zero),]
}

cat('Keep ', nrow(rna_seq), 'non-zero out of ', nrow(rna_seq_raw),'\n')
dim(rna_seq)
rna_seq$transcript_id = str_replace(rna_seq$gene, pattern = '[.][0-9]*', replacement = '')
row.names(rna_seq) = rna_seq$transcript_id

cat('Head of the rna_seq data:', '\n')
print(head(rna_seq[,1:10]))
#Find whether it's transcript or genes
if (length(grep('ENSG',rna_seq$gene[1])) > 0){
    duplicated_rows = duplicated(all.entrezgene$ensembl_gene_id)
    all.entrezgene = all.entrezgene[!duplicated_rows,]
    row.names(all.entrezgene) = all.entrezgene$ensembl_gene_id
}

#Change to the promoter regions. 
rna_seq[,c('chr','start','end', 'strand')] = all.entrezgene[rownames(rna_seq),c('chromosome_name','transcript_start', 'transcript_end', 'strand')]
rna_seq = rna_seq[complete.cases(rna_seq),]
tmp_data = rna_seq[,c('chr','start','end', 'strand')]
forward_strand = rna_seq$strand == 1
reverse_strand  = rna_seq$strand == -1
cat('Forward strand:', sum(forward_strand), 'Reverse ', sum(!forward_strand) ,'\n')

if(length(grep('genebody', batch_name)) > 0){#make it default to whole gene regions
    rna_seq[forward_strand,'start'] = tmp_data[forward_strand, 'start'] - 2000
    rna_seq[forward_strand,'end'] = tmp_data[forward_strand, 'end'] + 2000
    rna_seq[reverse_strand,'start'] = tmp_data[reverse_strand, 'start'] - 2000
    rna_seq[reverse_strand,'end'] = tmp_data[reverse_strand, 'end'] + 2000
}else{#Only the TSS regions
    rna_seq[forward_strand,'start'] = tmp_data[forward_strand, 'start'] - 5000
    rna_seq[forward_strand,'end'] = tmp_data[forward_strand, 'start'] + 5000
    rna_seq[reverse_strand,'start'] = tmp_data[reverse_strand, 'end'] - 5000
    rna_seq[reverse_strand,'end'] = tmp_data[reverse_strand, 'end'] + 5000
}



key_columns = c('chr','start','end', 'transcript_id', 'gene','strand')
rest_columns = setdiff(colnames(rna_seq), key_columns)
rna_seq_reordered = rna_seq[,c(key_columns, sample_cols)]
rna_seq_reordered$chr = paste0('chr', rna_seq$chr)

head(rna_seq_reordered)

##############Results#####################

cat("The number of transcript:", nrow(rna_seq_reordered), '\n')
print(table(rna_seq_reordered$chr))

library(dplyr)

chr_rows = grepl('chr[0-9X]+', rna_seq_reordered$chr)
cat('Missing transcripts', sum(!chr_rows), 'about', sum(!chr_rows)/nrow(rna_seq_reordered), '\n')

filtered_data = rna_seq_reordered[!chr_rows,]
head(filtered_data)
rna_seq_reordered2 = rna_seq_reordered[chr_rows,]
rna_seq_sorted = rna_seq_reordered2[with(rna_seq_reordered2, order(chr, start, end)), ]

head(rna_seq_sorted)

write.table(format(rna_seq_sorted, digits = 4), file = f_p('%s/rnaseq/transcript_data.bed', batch_output_dir), quote = F, sep = '\t', row.names = F)


#Write a small subset of know TFs.
tf_gene_id = read.table('./data/raw_data/rnaseq/tf_ensemble_id.txt', header = T, stringsAsFactors = FALSE)
rownames(tf_gene_id) = tf_gene_id$tf

tf_gene_expression = rna_seq_sorted[tf_gene_id$ensembl_gene_id,]

rownames(tf_gene_expression)

write.table(format(tf_gene_expression, digits = 4), f_p('./data/%s/rnaseq/tf_transcript_data.bed', batch_name), quote =FALSE, sep = '\t', row.names= FALSE)



####Create a random set of genes as control#######
set.seed(669)
faked_tfs =sample(rownames(rna_seq_sorted), size = nrow(tf_gene_expression))
faked_tf_expression = rna_seq_sorted[faked_tfs,]
head(faked_tf_expression)
rownames(faked_tf_expression) = rownames(tf_gene_expression)
write.table(format(faked_tf_expression, digits = 4), f_p('./data/%s/rnaseq/faked_tf_transcript_data.bed', batch_name), quote =FALSE, sep = '\t', row.names= FALSE)
dim(faked_tf_expression)

####Create a random set of miRNAs as control#####
miRNA_expression = read.table('./data/raw_data/miRNA/GD452.MirnaQuantCount.1.2N.50FN.samplename.resk10.txt', header = TRUE)
head10(miRNA_expression)
sample_cols = sort(grep('(NA|HG)[0-9]+', colnames(tf_gene_expression),value = T))
sample_cols = sort(intersect(sample_cols, colnames(miRNA_expression)))
random_mirna =sample(rownames(miRNA_expression), size = nrow(tf_gene_expression))
random_mirna_expression = miRNA_expression[random_mirna, sample_cols]




rownames(random_mirna_expression) = rownames(tf_gene_expression)
head10(random_mirna_expression)
dim(random_mirna_expression)

write.table(format(random_mirna_expression, digits = 4), f_p('./data/%s/rnaseq/random_miRNA_transcript_data.bed', batch_name), quote =FALSE, sep = '\t', row.names= FALSE)
cat('Writing Finished','\n')













