setwd('~/projects/expression_var/R/')
library(methods)
source('~/R/s_function.R', chdir = T)
library(futile.logger)
source('s_project_funcs.R')


GENE <- setRefClass("GENE",
                     fields = list( data = "data.frame",
                                   gene_name = "character",
                                   batch_name = 'character',
                                   chr_str = 'character'),
                     methods = list(

                         get_samples = function(){
                             sample_cols = sort(grep('(NA|HG)[0-9]+', colnames(data),value = T))
                             return (sample_cols)
                         },
                         
                         read_data = function(test_prefix=NULL) {
                             output_dir = f_p('./data/%s/', batch_name)
                             expression_file = f_p('%s/rnaseq/%s/%s.txt', output_dir, chr_str, gene_name)

                             
                             
                             data <<- read.table(expression_file, header = T, na.strings = 'NA')
                             
                             sample_cols = get_samples()
                             flog.info('Number of samples: %s', length(sample_cols))


                             if ('cor'  %in% colnames(data)){
                                 data <<- data[,c('chr', 'start', 'end', 'gene','feature',
                                                  'feature_start', 'feature_end', 'hic_fragment_id' ,
                                                  'type', 'cor', 'pair', sample_cols )]
                             }else{

                                 data <<- data[,c('chr', 'start', 'end', 'gene','feature',
                                                  'feature_start', 'feature_end', 'hic_fragment_id' ,
                                                  'type', sample_cols )]
                             }               
                             
                         },

                         change_expression =function(new_batch_name, batch_mode = 'TF') {
                             #Change the mRNA values to another normalization method,
                             #e.g. quantile to peer normalization.
                             if (!is.null(new_batch_name)){
                                 if (batch_mode != 'random'){
                                     new_data = read.table(f_p('./data/%s/rnaseq/GEUVADIS.Gene.DATA_MATRIX', new_batch_name), header = TRUE)
                                 }else{
                                     new_data = shuffle_gene_expression(new_batch_name)
                                 }
                                 rownames(new_data) = new_data$gene
                                 samples = get_samples()

                                 mRNA_index = which(data$type == 'gene')
                                 name_shared = intersect(samples , colnames(new_data))
                                 name_diff = setdiff(samples , colnames(new_data))
                                 flog.info('Set diff: %s', paste(name_diff, sep = '', collapse = '.'))
                                 data[mRNA_index, name_shared] <<- new_data[gene_name, name_shared]

          
                                 
                             }else{
                                 flog.info('New batch_name is null: %s', new_batch_name)
                             }
                         },

                         shuffle_gene_expression = function(loc_batch_name){
                             shuffle_file = f_p('./data/%s/rnaseq/GEUVADIS.GeneRandom.DATA_MATRIX', loc_batch_name)
                             if(!file.exists(shuffle_file)){
                                 loc_data = read.table(f_p('./data/%s/rnaseq/GEUVADIS.Gene.DATA_MATRIX', loc_batch_name), header = TRUE)
                                 permutate_names = sample(x =1:nrow(loc_data), size = nrow(loc_data))
                                 random_data = loc_data
                                 random_data$gene = loc_data[permutate_names, 'gene']
                                 write.table(random_data, file = shuffle_file, quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)
                                 flog.info('Write random expression file: %s', shuffle_file)
                             }else{
                                 random_data = read.table(shuffle_file, header = TRUE)
                             }
                             return (random_data)
                         },


                       subset_snps_in_tf_regions = function(batch_mode){

                           if ( !grepl( 'SNPinTF', batch_mode)) return (NULL)
                           
                           regulatory_regions = subset(data, type == 'promoter' | type == 'enhancer')[, c('chr', 'feature_start', 'feature_end', 'feature')]
                           tf_binding_regions = regulatory_regions[grep('H3K', regulatory_regions$feature, invert = T),]

                           SNP_regions = subset(data, type == 'SNP')[, c('chr', 'feature_start', 'feature_end', 'feature')]
                           
                           tf_bed = makeGRangesFromDataFrame(tf_binding_regions)
                           SNP_bed = makeGRangesFromDataFrame(SNP_regions)

                           matches = as.data.frame( findOverlaps(tf_bed, SNP_bed) )
                           removed_SNPs = setdiff(rownames(SNP_regions), rownames(SNP_regions[unique(matches[,2]),]))

                           if (batch_mode == 'randomSNPinTF'){
                               removed_SNPs = sample(rownames(SNP_regions), size = length(removed_SNPs))
                           }

                           flog.info('%s of %s SNPs are filtered out of the TF binding regions', length(removed_SNPs) , nrow(SNP_regions))
                           data <<- data[setdiff(rownames(data), removed_SNPs),]
                           
                       }
                         

                         
                         
                     ))

.DollarNames.GENE <- function(x, pattern){
    my_funcs = grep('callSuper|show|copy|export|field|getClass|getRefClass|import|trace|init|usingMethods',getRefClass(class(x))$methods(), value = T, invert = T) 
    grep(pattern, c( my_funcs, names(getRefClass(class(x))$fields())), value=TRUE)
}


library(RUnit)
options(run.main=FALSE)
if (getOption('run.main', default=TRUE)) {
    runTestFile('./test/t_gene_data_class.R',testFuncRegexp = '^t_.*')
}
