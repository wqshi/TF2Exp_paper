
#install.packages('bnlearn')
library(bnlearn)

expression_file = './data/rnaseq/gene_regulatory_region_feature_profile.chr1'
output_dir = './data/output/figures/'
expression_data = read.table(expression_file, header = TRUE)
head(expression_data)
cat('Number regulatory regions and empty fragment ids', '\n')
table(expression_data$type, expression_data$hic_fragment_id == '.')


table(unlist(table(expression_data$gene)))



#Get the expressed transcripts
sample_cols = grep('NA[0-9]+', colnames(expression_data),value = T)
rna_data = subset(expression_data, feature == 'RNASEQ')
rownames(rna_data) = rna_data$gene
expressed_transcripts = rownames(rna_data)[rowSums(rna_data[, sample_cols]) > 0]
cat('The number of expressed transcripts:', length(expressed_transcripts), '\n')



#Get the genes with at least one elements
gene_with_regulatory_elements = as.data.frame.matrix(table(expression_data$gene, expression_data$feature))
head(gene_with_regulatory_elements)
str(gene_with_regulatory_elements, max.level = 1)
regulated_genes_names = rownames(gene_with_regulatory_elements)[ rowSums(gene_with_regulatory_elements[,1:3]) > 5]
cat('The number of regulated transcripts:', length(regulated_genes_names), '\n')

#Get the inverstigeed genes
genes_names = intersect(expressed_transcripts, regulated_genes_names)
cat('The number of investigated transcripts:', length(genes_names), '\n')


length(genes_names)

library(caret)
for (i in 1:length(genes_names)){
    transcript_id = genes_names[i]
    transcript_data = expression_data[expression_data$gene == transcript_id,]



    train_data = transcript_data[,sample_cols]
    rownames(train_data) = make.names(paste0(transcript_data$type, '-',transcript_data$feature), unique = T)
    train_data_rmdup = train_data[!duplicated(train_data),]


##############################
#####train the BN model############
    train_data2 = as.data.frame(t(train_data_rmdup))
    correlated_features = findCorrelation(cor(train_data2), cutoff = 0.75)


    final_train_data  = train_data2[, -correlated_features]
    dim(final_train_data)
    if (ncol(final_train_data) < 50 & ncol(final_train_data) > 5){
        
        cat('Transcript', transcript_id, '-', i, 'has', ncol(final_train_data), 'features', '\n')
        res = fast.iamb(final_train_data)
        png(paste0(output_dir,transcript_data$gene[1], '.png'), width = 900, height = 900)
        plot(res, main = transcript_data$gene[1])
        fitted = bn.fit(res, final_train_data)
        dev.off()
        break
    }else{
        cat('Transcript', transcript_id, 'has', ncol(final_train_data), 'features', '\n')
    }


}


