
  
library(biomaRt)
#mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice"  )
#all.entrezgene <- unique( getBM(attributes = c("ensembl_gene_id", 'chromosome_name', 'strand', 'entrezgene'), values = "*", mart = mart) )
#dim(all.entrezgene)
#duplicated_rows = f_duplicated(all.entrezgene$entrezgene)
#write.table(all.entrezgene, file ='./data/raw_data/rnaseq/all.ensemble.genes.miRNA',sep='\t', quote = FALSE, row.names = FALSE )


all.entrezgene = read.table('./data/raw_data/rnaseq/all.ensemble.genes.miRNA',sep='\t', header = TRUE)




gene_miRNA_table = read.table('/homed/home/shi/expression_var/data/raw_data/miRNA/hsa_MTI.csv',sep=',', header = TRUE)
dim(gene_miRNA_table)

head(gene_miRNA_table)

colnames(gene_miRNA_table) = c('miRTarBase.id', 'miRNA_name', 'species', 'target_gene_name', 'target_entrez_id', 'gene_species', 'experiment', 'Experiment', 'Reference')

duplicated_rows = duplicated(x = gene_miRNA_table[,c('miRNA_name', 'target_entrez_id')])

subset_interaction_table = gene_miRNA_table[!duplicated_rows, c('miRNA_name', 'target_entrez_id')]
dim(subset_interaction_table)
head(subset_interaction_table)


ensemble_entrez_table = all.entrezgene[ complete.cases(all.entrezgene),c('ensembl_gene_id', 'entrezgene', 'chromosome_name')]
ensemble_entrez_table = ensemble_entrez_table[!duplicated(ensemble_entrez_table),]
dim(ensemble_entrez_table)
head(ensemble_entrez_table)

merged_df = merge(x = subset_interaction_table, y = ensemble_entrez_table, by.x = 'target_entrez_id', by.y = 'entrezgene' )

head(merged_df)

merged_df %>% filter(target_entrez_id == 1)
ensemble_entrez_table %>% filter(entrezgene == 1)
subset_interaction_table %>% filter(target_entrez_id == 1)

table(merged_df$chromosome_name)

length(unique(merged_df$miRNA_name))

write.table(merged_df, file ='/homed/home/shi/expression_var/data/raw_data/miRNA/miRNA_ensemble.txt', sep = '\t', row.names = F, quote = F)



