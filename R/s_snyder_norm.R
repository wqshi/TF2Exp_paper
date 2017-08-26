#The normalization method get from Michal Snyder.
source('~/R/s_function.R', chdir = TRUE)
library(matrixStats )
library(preprocessCore)
library(peer)
library(futile.logger)
source('s_project_funcs.R')
                                        #1) Variance stabilising transformation (i.e. asinh transformation) 
                                        # peaksMat is the count matrix for RNA transcripts (TPMs) obtained from Sailfish
                                        #batch_size = '462samples'

selected_samples = read.table('./data/raw_data/samples358.ped')
selected_samples$V1


#batch_size = '445samples'
batch_size = '358samples'
batch_name = f_p( '%s_sailfish', batch_size)
stablize_method = 'VarSqr' #Another way to stablized the variance of gene expression.

                                        #batch_name = '54samples'
batch_output_dir = f_p('./data/%s/', batch_name)

rna_seq_raw = read.csv(file = f_p('%s/rnaseq/GEUVADIS.Gene.DATA_MATRIX.sailfish', batch_output_dir) , sep=' ', header =TRUE)
rownames(rna_seq_raw) = rna_seq_raw$gene
rna_seq_raw$gene = NULL
if (grepl('358', batch_name)){
    rna_seq_raw = rna_seq_raw[, colnames(rna_seq_raw) %in% selected_samples$V1]
}


rna_seq_raw_rmZero = f_filter_NA_rm_zeroVar(rna_seq_raw)

gene = data.frame( gene = row.names(rna_seq_raw_rmZero) )
write_data = cbind(gene, rna_seq_raw_rmZero)


write.table( c(genes=nrow(rna_seq_raw_rmZero), samples = ncol(rna_seq_raw_rmZero)   ), './result/stats/number_of_gene.txt', quote = F)
write.table(write_data, file = f_p('%s/rnaseq/GEUVADIS.Gene.DATA_MATRIX', batch_output_dir), quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)


raw_index=f_check_na_rows(rna_seq_raw_rmZero)

peaksMat = as.matrix(rna_seq_raw_rmZero)

peaksMat_asinh = asinh(peaksMat) 
f_check_na_rows(peaksMat_asinh)

library(ggplot2)


p<- qplot(peaksMat[,1], (peaksMat_asinh[,1])^2) + xlab('Raw gene expression (TPM)') + ylab('Varaince stabilized (asinh(x))')
ggsave('./data/r_results/expression_stabilization.tiff', plot = p)

pp<-qplot(peaksMat[,1], (log2(peaksMat[,1])) ) + xlab('Raw gene expression (TPM)') + ylab('Varaince stabilized (log + square)')
ggsave('./data/r_results/expression_stabilization_log.tiff', plot = pp)

hist((peaksMat_asinh[,1])^2)
hist((peaksMat_asinh[,1]))
                                        # 2) Standardisation:

rMs = rowMeans(peaksMat_asinh, na.rm=TRUE)
rSds = rowSds(peaksMat_asinh, na.rm=TRUE)
peaksMat_asinh_std = (peaksMat_asinh - rMs)/rSds
f_check_na_rows(peaksMat_asinh_std)
head10(peaksMat_asinh_std)

target_gene = 'ENSG00000063515.2'
#peaksMat_asinh_std[target_gene,1:10]

                                        # 3) remove genes with zero SD
peaksMat_asinh_std = peaksMat_asinh_std[rSds>0,]
peaksMat_asinh_std = peaksMat_asinh_std[ rowSums(is.na(peaksMat_asinh_std))  == 0, ]
f_check_na_rows(peaksMat_asinh_std)

                                        # 4) quantile normalize samples
peaksMat_asinh_std_qn = normalize.quantiles(peaksMat_asinh_std)
rownames(peaksMat_asinh_std_qn) = rownames(peaksMat_asinh_std)
colnames(peaksMat_asinh_std_qn) = colnames(peaksMat_asinh_std)


                                        #peaksMat_asinh_std_qn[selected_genes, selected_indivs]
                                        #peaksMat_asinh_std[selected_genes, selected_indivs]
gene = data.frame( gene = row.names(peaksMat_asinh_std_qn) )
write_data1 = cbind(gene, peaksMat_asinh_std_qn)
                                        #write.table(write_data1, file = './data/462samples_sailfish_quantile/rnaseq/GEUVADIS.Gene.DATA_MATRIX', quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)




##Read covariate
sample_info = read.table(f_p('./data/462samples/chr_vcf_files/integrated_call_samples_v3.20130502.ALL.panel'), header = TRUE, row.names = 1)
covariate=sample_info[colnames(peaksMat_asinh_std_qn), c('pop', 'gender')]



###########Approach remove all the covariate######

snyder_original_dir = f_p('./data/%s_snyder_original/rnaseq/', batch_size)
dir.create(snyder_original_dir, recursive = T)
na_index=f_check_na_rows(peaksMat_asinh_std_qn)
flog.info('Before peer, NA rows %s', length(na_index))


                                        #Regressively remove the batch effects.
f_cap_matrix <- function(loc_data, low_threshold, up_threshold){
    cap_data =loc_data
    cap_data[cap_data < low_threshold] = low_threshold
    cap_data[cap_data > up_threshold] = up_threshold
    return (cap_data)
}


set.seed(191)
#f_heatmap_genes(f_cap_matrix(peaksMat_asinh_std_qn, -5, 5), f_p('%s/gene_exp_heatmap_peer%s.tiff', snyder_original_dir, 0), seed_number = 11)



peer_range = c(10)

for (peerFactor in peer_range){
    model = PEER()
    PEER_setPhenoMean(model,t(as.matrix(peaksMat_asinh_std_qn)))
    PEER_setAdd_mean(model, FALSE)
    PEER_setNk(model,peerFactor)
    PEER_setNmax_iterations(model, 1000)
    PEER_update(model)

    tiff(filename = f_p('%s/peer_no_cov%s.tif', snyder_original_dir, peerFactor))
    #PEER_plotModel(model)
    alpha = PEER_getAlpha(model)
    plot(1/alpha,xlab="Factors",ylab="Variance of factor weights", type="b", col="blue", lwd=4, xaxp=c(1,length(alpha), length(alpha)-1),
         cex.lab=1.5, cex.axis = 1.5)
    
    dev.off()

    f_peer_variance_plot(model)
    
    factors = PEER_getX(model)
    weights = PEER_getW(model)
    precision = PEER_getAlpha(model)
    residuals = t(PEER_getResiduals(model))
    peaksMat_peer = (residuals-rowMeans(residuals))/rowSds(residuals)
    peaksMat_norm = normalize.quantiles(peaksMat_peer)
    colnames(peaksMat_norm)=colnames(peaksMat_asinh_std_qn)
    rownames(peaksMat_norm)=rownames(peaksMat_asinh_std_qn)            

    gene = data.frame( gene = row.names(peaksMat_norm) )
    write_data = cbind(gene, peaksMat_norm)
    head(write_data[,1:10])
    #f_heatmap_genes(f_cap_matrix(peaksMat_norm, -5, 5), f_p('%s/gene_exp_heatmap_peer%s.tiff', snyder_original_dir, peerFactor), seed_number=11)
}
flog.info(snyder_original_dir)
write.table(write_data, file = f_p('%s/GEUVADIS.Gene.DATA_MATRIX', snyder_original_dir), quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)
head(write_data[1:10, 1:4])


#stop() #For testing the plot_peer function



########Explicitly remove the covariate######
snyder_norm_dir = f_p('./data/%s_snyder_norm/rnaseq/', batch_size)
dir.create(snyder_norm_dir, recursive = T)

peer_range = seq(from = 40, to =60, by = 10)
peer_range = c(40)

for (peerFactor in peer_range){
    flog.info('Hidden factor %s', peerFactor)
    model = PEER()
    na.index=f_check_na_rows(peaksMat_asinh_std_qn)

    dim(peaksMat_asinh_std_qn)

    colnames(peaksMat_asinh_std_qn) == rownames(covariate)

    PEER_setNk(model,peerFactor)
    PEER_setPhenoMean(model,t(as.matrix(peaksMat_asinh_std_qn)))

    ##Peer needs numeric values for covariates.
    source('s_gene_regression_fun.R')
    head(covariate)
    table(covariate$pop)
    covariate$pop = as.character(covariate$pop)

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
    PEER_setNmax_iterations(model, 1000)
    PEER_update(model)


    tiff(filename = f_p('%s/peer%s.tif', snyder_norm_dir, peerFactor), width = 1080, height = 600)
    #PEER_plotModel(model)
    alpha = PEER_getAlpha(model)
    plot(1/alpha,xlab="Factors",ylab="Variance of factor weights", type="b", col="blue", lwd=4, xaxp=c(1,length(alpha), length(alpha)-1),
         cex.lab=1.5, cex.axis = 1.5)
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

    gene = data.frame( gene = row.names(peaksMat_norm) )
    cov_write_data = cbind(gene, cov_peer_norm)
    head(write_data[,1:10])
    #f_heatmap_genes(cov_peer_norm, f_p('%s/gene_exp_heatmap.tiff', snyder_norm_dir))
    library(futile.logger)
    flog.info(snyder_norm_dir)
    write.table(cov_write_data, file = f_p('%s/GEUVADIS.Gene.DATA_MATRIX.peer%s', snyder_norm_dir, peerFactor), quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)
}



