source('s_gene_data_class.R')



CEU60_samples_df = read.table('./data/raw_data/CEU60/E-MTAB-197.sdrf.txt', sep ='\t', header = TRUE)
CEU60_samples=unique(CEU60_samples_df$Factor.Value.INDIVIDUAL.)


YRI69_samples_df = read.table('./data/raw_data/YRI69/E-GEOD-19480.sdrf.txt', sep ='\t', header = TRUE)
YRI69_samples=unique(YRI69_samples_df$FactorValue..CELL.LINE.)



batch_name = '462samples_log_quantile'
gene = 'ENSG00000269103.1'
chr_str = 'chr22'
test_gene <- GENE(data = data.frame(), gene_name = gene, chr_str = chr_str, batch_name = batch_name)
test_gene$read_data()
samples370 = test_gene$get_samples()

unique_CEU_samples = setdiff(CEU60_samples, samples370)
unique_YRI_samples = setdiff(YRI69_samples, samples370)


###GEUVIS samples.
geuvis_samples = read.table('/homed/home/shi/expression_var/data/462samples/output/sample.list')
colnames(geuvis_samples) = c('id', 'sample')
head(geuvis_samples)


###Why many papers say it's 420 overlap with 1000 Genome project###
#' #1. The number individuals in other RNAseq data.


kg_file = '/homed/home/shi/expression_var/python/data/raw_data/wgs/1kg/integrated_call_samples_v3.20130502.ALL.panel'
kg_samples = read.table(kg_file, header = T)

cat(nrow(kg_samples), '1kg samples \n')

cat('1kg Overlaped', length(intersect(geuvis_samples$sample, kg_samples$sample)), 'with geuvis \n')

overlapped_samples=intersect(geuvis_samples$sample, kg_samples$sample)
write.table(overlapped_samples, './data/raw_data/samples445.ped', quote=F, row.names =F, col.names = F)



non_YRI_samples=intersect(geuvis_samples$sample, subset(kg_samples, pop != 'YRI')$sample)
length(non_YRI_samples)
write.table(non_YRI_samples, './data/raw_data/samples358.ped', quote=F, row.names =F, col.names = F)




unique_CEU_samples = intersect(kg_samples$sample, setdiff(CEU60_samples, geuvis_samples$sample))
cat(length(unique_CEU_samples), 'of CEU60 can be used for testing \n')
unique_YRI_samples2 = intersect(kg_samples$sample, setdiff(YRI69_samples, geuvis_samples$sample))
cat(length(unique_YRI_samples2), 'of YRI can be used for testing \n')


#' #Check the number of download files:
loire.part = read.table('/homed/home/shi/expression_var/sailfish/data/fastq250/loire.part')
clustdell.part = read.table('/homed/home/shi/expression_var/sailfish/data/fastq250/clustdell.part')
combined_df = rbind(loire.part, clustdell.part)
download_data = unique(grep('fastq.gz', combined_df$V1, value = T))
write.table(download_data, '/homed/home/shi/expression_var/sailfish/data/finish.20160710.part', quote = F, col.names = F)


#' #727 individuals from the microarray data
array800 = read.table('/homed/home/shi/expression_var/data/800samples/rnaseq/727samples.list')
cat(length(intersect(kg_samples$sample, array800$V1)), 'out of 727 samples are overlapped with 1kg project', '\n')
cat(length(intersect(array800$V1, geuvis_samples$sample)), 'out of 727 samples are overlapped with geuvis', '\n')
cat(length( intersect(kg_samples$sample, setdiff(array800$V1, geuvis_samples$sample))), 'out of 727 samples are overlapped with kg but not with geuvis', '\n')

test_samples=intersect(kg_samples$sample, setdiff(array800$V1, geuvis_samples$sample))
rownames(kg_samples) = kg_samples$sample
kg_samples$pop = as.character(kg_samples$pop)
#' #distribution of available test samples.
print(table(kg_samples[test_samples,'pop']))


write.table(intersect(kg_samples$sample, array800$V1), './data/800samples/output/samplesOverlap1KG.list', quote = F, col.names = F)


write.table(unique(c(intersect(kg_samples$sample, array800$V1), intersect(geuvis_samples$sample, kg_samples$sample)), 'NA12878'), './data/800samples/output/AllSamplesOverlap1KG.list', quote = F, col.names = F)

###Check the SNP difference between two versions.
snp_358 = read.table('/homed/home/shi/expression_var/data/raw_data/wgs/1kg/additive_358samples/chr22.vcf.snps')
snp_445 = read.table('/homed/home/shi/expression_var/data/raw_data/wgs/1kg/additive_445samples/chr22.vcf.snps')

head(snp_358)

length( setdiff(snp_445$V1, snp_358$V1))
length( setdiff(snp_358$V1, snp_445$V1))


#########
#Parse additonal test data.








