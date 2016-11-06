
source('~/R/s_function.R', chdir = TRUE)
source('s_project_funcs.R')

f_read_expression_data <- function(input_file){
    exp_data = read.table(input_file, sep = ' ', header = T)
    rownames(exp_data) = str_replace(exp_data$gene, '[.].*', '')
    exp_data$gene = NULL
    return (exp_data)
}

exp_800 = f_read_expression_data('./data/800samples/rnaseq/GEUVADIS.Gene.DATA_MATRIX')
exp_445 = f_read_expression_data('./data/445samples_sailfish/rnaseq/GEUVADIS.Gene.DATA_MATRIX')


dim(exp_445)
head10(exp_445)
colnames(exp_800)


shared_indvs = intersect( colnames(exp_800), colnames(exp_445) )
shared_genes = intersect( rownames(exp_800), rownames(exp_445)  )


loc_gene = shared_genes[1]
sample_info = read.table(f_p('./data/462samples/chr_vcf_files/integrated_call_samples_v3.20130502.ALL.panel'), header = TRUE, row.names = 1)
head(sample_info)

#shared_genes = shared_genes[1:10]

test_samples = sample_info[shared_indvs,]
loc_pop = 'CEU'
loc_gene = shared_genes[1]
for(loc_pop in unique(test_samples$pop)){
    
    pop_shared_indiv = rownames(subset(test_samples, pop == loc_pop))
    cat('========', loc_pop, ':', length(pop_shared_indiv), ' individuals ========')
    result_list = rep(0, length(shared_genes))
    names(result_list) = shared_genes
    for (loc_gene in shared_genes){
        #cat(loc_gene,'\n')
        result_list[loc_gene]=cor(  x = as.numeric(exp_445[loc_gene, pop_shared_indiv]), y = as.numeric( exp_800[loc_gene, pop_shared_indiv] ), method = 'spearman', use = 'complete.obs')    
    }

    hist(result_list[shared_genes])
    cat('Total', length(shared_genes),'Mean correlation', mean(result_list, na.rm = TRUE),'\n')

}


individuals_results <- rep(0, times = length(shared_indvs))
names(individuals_results) = shared_indvs

for ( loc_indiv in shared_indvs ){
    individuals_results[loc_indiv]=cor(exp_445[shared_genes, loc_indiv], exp_800[shared_genes, loc_indiv],  method = 'spearman', use = 'complete.obs' )
}

cat('Total', length(shared_indvs),'Mean correlation', mean(individuals_results, na.rm = TRUE),'\n')

