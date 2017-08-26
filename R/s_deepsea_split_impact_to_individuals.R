source('~/R/s_function.R', chdir = TRUE)
library(stringr)
source('s_project_funcs.R')
library(foreach)
library(doMC)
registerDoMC(cores=4)
source('s_gene_regression_fun.R')
#Input
#batch_name = '462samples_sailfish'
#batch_name = '800samples'
#batch_name = '445samples_diff'
batch_name = '445samples_region'
batch_name = '800samples_region'

###################
##The rest code##

chr_str = 'chr22'


f_split_one_chromosom <- function(chr_str, deepsea_dir, batch_name, test_flag = FALSE){
         
    cat('============', chr_str, '=========\n')
    result_dir = f_p('/homed/home/shi/expression_var/data/%s/deep_result/all/%s/evalue/', batch_name, chr_str)
    vcf_dir = f_p('/homed/home/shi/expression_var/data/%s/chr_vcf_files/%s/', batch_name, chr_str)
    dir.create(result_dir, recursive = T)
        
    chr_impact = read.table(f_p('%s/%s.diff.gz', deepsea_dir, chr_str), sep = ',', stringsAsFactors=FALSE)
        
    colnames(chr_impact) = chr_impact[1,]
    chr_impact = chr_impact[chr_impact$chr != 'chr',]
    rownames(chr_impact) = paste(chr_impact$chr, chr_impact$pos, chr_impact$ref, chr_impact$alt, sep = ':')
    head10(chr_impact)

    col_names = str_replace_all(colnames(chr_impact),'None[.][0-9]+$', 'None')

    colnames(chr_impact)

    cat('Number of variations in total:', nrow(chr_impact), '\n')
    head10(chr_impact)

    individual_vcfs_list = grep('vcf.gz', list.files(vcf_dir), value = TRUE)
    cat(length(individual_vcfs_list), 'individuals', '\n')
    individual_vcf = individual_vcfs_list[1]
    if (test_flag == TRUE){
        individual_vcfs_list=individual_vcfs_list[1:2]
    }
    individual_vcf = individual_vcfs_list[2]
    
    for (individual_vcf in individual_vcfs_list){
        
        individual_id = str_replace(individual_vcf, '.vcf.gz', '')
        cat(individual_id, '\n')
        
        result <- try(
                        individual_genotype <- read.table(f_p('%s/%s', vcf_dir, individual_vcf), sep = '\t'),
                        silent=TRUE
                       )
    
        if (class(result) == "try-error"){
            print(result)
            cat('Error in ', individual_id, '\n')
            next
        }

        individual_genotype = individual_genotype[!duplicated(individual_genotype),]
       
        colnames(individual_genotype) = c('chr', 'pos', 'genotype', 'ref', 'alt')
        individual_genotype = individual_genotype[grepl('^[ATGCatgc]+$', individual_genotype$alt),] 
        
        rownames(individual_genotype) = with(individual_genotype, paste(chr, pos, ref, alt, sep = ':'))
        
        individual_impact = chr_impact[rownames(individual_genotype),]
        individual_impact$name = individual_genotype[rownames(individual_impact), 'genotype']
        f_ASSERT(sum(individual_impact$name == individual_genotype$genotype)==nrow(individual_genotype), 'Genotype mis-match')
        
        #head10(individual_impact)

        gz1 <- gzfile(f_p("%s/%s.diff.gz", result_dir, individual_id ), "w")
        write_data = rbind(col_names, individual_impact)
        flog.info('Remove %s NA rows out of %s', nrow(write_data) - sum(complete.cases(write_data)), nrow(write_data))
        write.table(write_data[complete.cases(write_data),], gz1, sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE)
        close(gz1)
    }
}


#chr_list = c('chr22', 'chr10', 'chr15')
chr_list = paste0('chr', c(1:22, 'X'))

for (batch_name in c('445samples_region')){

    if (batch_name == '800samples_region'){
        deepsea_dir = '/homed/home/shi/expression_var/data/800samples_region/deep_result/all/chr800/evalue/'
    }else{
        deepsea_dir = f_p('/homed/home/shi/expression_var/data/%s/deep_result/all/chrMerge2/evalue/', batch_name)
    }


    for (chr_str in chr_list){
        f_split_one_chromosom(chr_str, deepsea_dir, batch_name, test_flag = FALSE)
        data_dir=f_p('/home/shi/expression_var/data/%s/deep_result/all/%s/evalue/', batch_name, chr_str)
        system(f_p('scp /homed/%s/*.gz $CLUST:%s', data_dir, data_dir))
    }

}

#foreach( chr_str = paste0('chr', c(22:1,'X')) ) %dopar% f_split_one_chromosom(chr_str, deepsea_dir, batch_name)

#foreach( chr_str = paste0('chr', 21:1), .combine=c ) %do% f_split_one_chromosom(chr_str, deepsea_dir)
#print(results) 






