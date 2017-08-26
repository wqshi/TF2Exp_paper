library(peer)
source('~/R/s_function.R', chdir = TRUE)
#install.packages('preprocessCore')
#source("http://bioconductor.org/biocLite.R")
#biocLite("preprocessCore")
source('s_project_funcs.R')
library(preprocessCore)
library(matrixStats)
library(reshape2)
source('s_gene_regression_fun.R')


##Read the ChIP-seq data

input_file = '/homed/home/shi/expression_var/R/data/raw_data/tf/pu1_47indivs.txt'
factor_name = 'pu1'

input_file = '/homed/home/shi/expression_var/R/data/raw_data/tf/pu1_47indivs.txt'
factor_name = 'pu1'


f_preprocess_hommer_output <- function(input_file, factor_name){

    binding_data_raw = read.table(input_file, sep = '\t', fill = TRUE,quote = "")
    binding_data = binding_data_raw
    head10(binding_data_raw)
    col_names  =gsub('.*(NA[0-9].*)[/].*','\\1',binding_data[1,])
    col_names[1] = 'PeakID'
    table(binding_data$Strand)
    colnames(binding_data) = col_names
    binding_data = subset(binding_data, Strand == '+')


    binding_data$PeakID = paste0( binding_data$Chr, '_', as.numeric(binding_data$Start)-1, '_', binding_data$End )

    sample_cols = grep('(NA|HG)[0-9]+', colnames(binding_data),value = T)
    binding_matrix = binding_data[,c('PeakID', sample_cols)]
    head10(binding_matrix)
    write.table(binding_matrix, file = f_p('/homed/home/shi/expression_var/R/data/raw_data/CEU47/%s_raw.txt', factor_name), quote =F, sep = '\t')
    return (binding_matrix)


}


f_get_mean_of_replicates <- function(input_data){

    
    sample_cols = grep('NA', colnames(input_data), value = T)
    sample_ids = unique(str_replace(sample_cols, '_Rep.*', ''))

    new_data = data.frame(matrix(0, nrow = nrow(input_data), ncol = length(sample_ids)+1))
    colnames(new_data) = c('PeakID', sample_ids)
    new_data$PeakID = input_data$PeakID
    loc_sample = sample_ids[1]

    for(loc_sample in sample_ids){
        match_cols = grep(loc_sample, sample_cols, value = T)
        cat(match_cols,'\n')
        if (length(match_cols) > 1){
            new_data[,loc_sample] = apply(data.matrix(input_data[,match_cols]), 1, mean)
        }else{
            new_data[,loc_sample] = data.matrix(input_data[,match_cols])
        }
    }

    return (new_data)
}


f_peer_normalize_data <- function(input_data, id_col, factor_name, output_dir, peer_factor = NULL, stablize_exp = FALSE, debug = TRUE){

    if (debug == TRUE){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##

    }

    rownames(input_data) = input_data[[id_col]]

    input_data[[id_col]] = NULL

    rna_seq_raw_rmZero = f_filter_NA_rm_zeroVar(input_data)

    dim(rna_seq_raw_rmZero)

    peaksMat = data.matrix(rna_seq_raw_rmZero)

    head(peaksMat)

    if(stablize_exp == TRUE){
        peaksMat_asinh = asinh(peaksMat) 
    }
    
    peaksMat_asinh = peaksMat


    rMs = rowMeans(peaksMat_asinh, na.rm=TRUE)
    rSds = rowSds(peaksMat_asinh, na.rm=TRUE)
    peaksMat_asinh_std = (peaksMat_asinh - rMs)/rSds
    f_check_na_rows(peaksMat_asinh_std)
    head10(peaksMat_asinh_std)

                                        #peaksMat_asinh_std[target_gene,1:10]

                                        # 3) remove genes with zero SD
    peaksMat_asinh_std = peaksMat_asinh_std[rSds>0,]
    peaksMat_asinh_std = peaksMat_asinh_std[ rowSums(is.na(peaksMat_asinh_std))  == 0, ]
    f_check_na_rows(peaksMat_asinh_std)

                                        # 4) quantile normalize samples
    peaksMat_asinh_std_qn = normalize.quantiles(peaksMat_asinh_std)
    rownames(peaksMat_asinh_std_qn) = rownames(peaksMat_asinh_std)
    colnames(peaksMat_asinh_std_qn) = colnames(peaksMat_asinh_std)




    sample_info = read.table(f_p('./data/462samples/chr_vcf_files/integrated_call_samples_v3.20130502.ALL.panel'), header = TRUE, row.names = 1)
    shared_indivs = intersect( rownames(sample_info), colnames(rna_seq_raw_rmZero) )
    peaksMat_asinh_std_qn = peaksMat_asinh_std_qn[, shared_indivs]
    covariate=sample_info[colnames(peaksMat_asinh_std_qn), c('pop', 'gender')]
    

    if (is.null(peer_factor)){
        peer_range = seq(from = 10, to =15, by = 5)
    }else{
        peer_range = c(peer_factor)
    }
    #peerFactor = 40

    for (peerFactor in peer_range){
        flog.info('Hidden factor %s', peerFactor)
        model = PEER()
        na.index=f_check_na_rows(peaksMat_asinh_std_qn)

        dim(peaksMat_asinh_std_qn)

        f_ASSERT( all(colnames(peaksMat_asinh_std_qn) == rownames(covariate)), "Rownames don't match")

        PEER_setNk(model,peerFactor)
        PEER_setPhenoMean(model,t(as.matrix(peaksMat_asinh_std_qn)))

        ##Peer needs numeric values for covariates.

        head(covariate)
        table(covariate$pop)
        #covariate$pop = as.character(covariate$pop)

        if(length(unique(covariate$pop)) == 1){
            covariate$pop = NULL
        }

        
        head(covariate)
        dummies <- dummyVars(as.formula(' ~ .'), data = covariate)
        
        cov_dummy=as.data.frame(predict(dummies, newdata = covariate))

                
        print(colnames(cov_dummy))
        ##cov_dummy$gender.male = NULL
        cov_dummy$genderfemale = NULL
        ##cov_dummy$popCEU = NULL
        colSums(cov_dummy)
        ##PEER_setCovariates(model, as.matrix(cov_dummy[,c('gender.female', 'gender.male')]))
        PEER_setCovariates(model, as.matrix(cov_dummy))
        PEER_getCovariates(model)
        ##PEER_setAdd_mean(model, TRUE)
        PEER_setNmax_iterations(model, 300)
        converge_state = PEER_update(model)
        flog.info('Converge state', converge_state)

        tiff(filename = f_p('%s/peer%s_%s.tif', output_dir, peerFactor, factor_name), width = 1080, height = 1080)
        PEER_plotModel(model)
        dev.off()

        cov_factors = PEER_getX(model)
        head(cov_factors)

        cov_weights = PEER_getW(model)
        head(cov_weights)
        cov_precision = PEER_getAlpha(model)
        cov_residuals = t(PEER_getResiduals(model))
        cov_peer = (cov_residuals-rowMeans(cov_residuals))/rowSds(cov_residuals)
        cov_peer_norm = normalize.quantiles(cov_peer)
        colnames(cov_peer_norm)=colnames(peaksMat_asinh_std_qn)
        rownames(cov_peer_norm)=rownames(peaksMat_asinh_std_qn)            

        gene = data.frame( gene = row.names(peaksMat_asinh_std_qn) )
        colnames(gene) = id_col
        cov_write_data = cbind(gene, cov_peer_norm)
        
        ##f_heatmap_genes(cov_peer_norm, f_p('%s/gene_exp_heatmap.tiff', snyder_norm_dir))
        library(futile.logger)
        flog.info(output_dir)
        write.table(cov_write_data, file = f_p('%s/peer%s_%s.txt', output_dir, peerFactor, factor_name), quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)
    }

}




id_col = 'PeakID'
output_dir = '/homed/home/shi/expression_var/data/raw_data/CEU47/'

###PU1
#factor_name = 'PU1'
#input_data = f_preprocess_hommer_output(input_file, factor_name)
#head10(input_data)


factor_name = 'CTCF'
input_data = f_preprocess_hommer_output('/homed/home/shi/expression_var/R/data/raw_data/CEU_CTCF/output.file', factor_name)
mean_input = f_get_mean_of_replicates(input_data)






#Test


hist(as.numeric(input_data[,2]))
#f_peer_normalize_data(input_data, id_col, factor_name, output_dir, debug = F)
#f_peer_normalize_data(input_data, id_col, factor_name, output_dir, peer_factor = 9, debug = F)


dim(mean_input)


##CTCF
f_peer_normalize_data(mean_input, id_col, factor_name, output_dir, peer_factor = NULL, debug = F)
f_peer_normalize_data(mean_input, id_col, factor_name, output_dir, peer_factor = 9, debug = F)




stop()
##Normalize expression data
expression_data = expression_47indivs = read.table('/homed/home/shi/expression_var/data/raw_data/CEU47/GEUVADIS.Gene.DATA_MATRIX', header = T, sep = ' ')

head(expression_data)
#f_peer_normalize_data(expression_data, id_col = 'gene', factor_name='gene', output_dir, peer_factor = NULL, stablize_exp = T, debug = F)
f_peer_normalize_data(expression_data, id_col = 'gene', factor_name='gene', output_dir, peer_factor = 10, stablize_exp = T, debug = F)











