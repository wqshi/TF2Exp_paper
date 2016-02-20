library(peer)
source('~/R/s_function.R', chdir = TRUE)
#install.packages('preprocessCore')
#source("http://bioconductor.org/biocLite.R")
#biocLite("preprocessCore")
library(preprocessCore)
library(matrixStats)
library(reshape2)

f_row_scale <- function(df){
    return (t(scale(t(df))))
}

f_quantile <- function(input_data){

    data_quantile = normalize.quantiles(as.matrix(input_data))
    colnames(data_quantile) = colnames(input_data)
    rownames(data_quantile) = rownames(input_data)
    return (data_quantile)
}



batch_name = '462samples_sailfish'
#batch_name = '54samples_evalue'
batch_output_dir = f_p('./data/%s/', batch_name)

rna_seq_raw = read.csv(file = f_p('%s/rnaseq/GEUVADIS.Gene.DATA_MATRIX', batch_output_dir) , sep=' ', header =TRUE)
rownames(rna_seq_raw) = rna_seq_raw$gene
rna_seq_raw$gene = NULL


selected_genes = c('ENSG00000177663.9', 'ENSG00000235478.1')
selected_indivs = c('HG00145', 'HG00310', 'HG00171')



########Check the data ####################
#rna_seq_raw[is.na(rna_seq_raw)] = 0 #This is wrong
sample_cols = grep('(NA|HG)[0-9]+', colnames(rna_seq_raw),value = T)
head10(rna_seq_raw)

set.seed(11)
#Check the assumption for quantile normalization
random_samples = sample(colnames(rna_seq_raw), size = 15)
long_exp = melt(rna_seq_raw[, random_samples])
head(long_exp)
sample_plot <- ggplot(long_exp, aes(x = value, group = variable, color = variable)) + geom_density() + xlim(c(0,50))





##Case1: Quantile normalization the raw gene expression data.
rna_seq_quantile = normalize.quantiles(as.matrix(rna_seq_raw))
colnames(rna_seq_quantile) = colnames(rna_seq_raw)

cat('NA rows in the data:',  sum(rowSums(is.na(rna_seq_quantile)) > 0), '\n')

head(rna_seq_quantile[,1:10])
head(rna_seq_raw[,1:10])

rna_seq = scale(rna_seq_raw)

rna_seq_scale = f_row_scale(rna_seq_quantile)

head(rna_seq_scale[,1:10])
head(rna_seq[,1:10])
dim(rna_seq_scale)


gene = data.frame( gene = rownames(rna_seq_raw) )
write_data = cbind(gene, rna_seq_scale)
rownames(write_data) = write_data$gene
#head(write_data[,1:10])

cat('Check it right:', '\n')
#ENSG00000177663.9  0.2461960 -0.6264528 -0.7476451
#ENSG00000235478.1 -0.3305492  4.4925357 -0.3305492
print(write_data[selected_genes, selected_indivs])
write.table(write_data, file = './data/462samples_quantile/rnaseq/GEUVADIS.Gene.DATA_MATRIX', quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)



#################################
#################################
#Case 2: quantile_rmNA: remove the gene with NA expression and zero expression.


rnaseq_rmNA = rna_seq_raw[rowSums(is.na(rna_seq_raw)) == 0,]
cat('Keep', nrow(rnaseq_rmNA), 'out of ', nrow(rna_seq_raw) , 'none-NA genes \n')

library(stringr)
sex.entrezgene = read.table('./data/raw_data/rnaseq/sex.ensemble.genes',sep='\t', header = TRUE)
row.names(sex.entrezgene) = sex.entrezgene$ensembl_transcript_id

rnaseq_rmNA$gene = str_replace( row.names(rnaseq_rmNA) )

grep(sex_pattern, rownames(rnaseq_rmNA), invert = TRUE)
head10(rnaseq_rmNA)



rnaseq_clean = rnaseq_rmNA[(rowSds(as.matrix(rnaseq_rmNA)) != 0),]
cat('Keep', nrow(rnaseq_clean), 'out of ', nrow(rnaseq_rmNA) , 'zero expressed genes \n')







head10(rnaseq_clean)

clean_quantile = f_quantile(rnaseq_clean)
clean_quantile_scale = f_row_scale(clean_quantile)
f_check_na(clean_quantile_scale)
head10(clean_quantile_scale)
gene = data.frame( gene = rownames(clean_quantile_scale) )
write_data = cbind(gene, clean_quantile_scale)
rownames(write_data) = write_data$gene
#head(write_data[,1:10])
write.table(write_data, file = './data/462samples_quantile_rmNA/rnaseq/GEUVADIS.Gene.DATA_MATRIX', quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)




################################
################################
#Case 3: variance stablize before quantile_rmNA

peaksMat = as.matrix(rnaseq_clean)
peaksMat_asinh = asinh(peaksMat) 

head10(peaksMat_asinh)
asinh_quantile = f_quantile(peaksMat_asinh)
asinh_scale = f_row_scale(asinh_quantile)

gene = data.frame( gene = rownames(asinh_scale) )
asinh_write_data = cbind(gene, asinh_scale)
write.table(asinh_write_data, file = './data/462samples_var_quantile/rnaseq/GEUVADIS.Gene.DATA_MATRIX', quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)



stop()
#################################
#################################




PEER_plotModel <- function(model){
    par(mfrow=c(2,1))
    bounds = PEER_getBounds(model)
    vars = PEER_getResidualVars(model)
    par(mar=c(5,4,4,5)+.1)
    plot(bounds, type="l", col="red", lwd=2, xlab="Iterations", ylab="Lower bound")
    par(new=TRUE)
    plot(vars,,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
    axis(4)
    mtext("Residual variance",side=4,line=3)
    legend("right",col=c("red","blue"),lty=1,legend=c("Lower bound","Residual variance"))
    alpha = PEER_getAlpha(model)
    plot(alpha,xlab="Factors",ylab="Inverse variance of factor weights", type="b", col="blue", lwd=4, xaxp=c(1,length(alpha), length(alpha)-1))
}



peaksMat = as.matrix(rna_seq_raw)
peaksMat_asinh = asinh(peaksMat) 

peer_model = PEER()



PEER_setPhenoMean(peer_model, t(peaksMat_asinh))
dim(PEER_getPhenoMean(peer_model))
PEER_setNk(peer_model,30)
PEER_setNmax_iterations(peer_model, 1000)

PEER_update(peer_model)

tiff(filename = f_p('./data/462samples_peer/rnaseq/peer.tif'))
PEER_plotModel(peer_model)
dev.off()


model = peer_model
factors = PEER_getX(model)


#head(rna_seq_raw)


dim(factors)
weights = PEER_getW(model)
dim(weights)
precision = PEER_getAlpha(model)
dim(precision)
residuals = PEER_getResiduals(model)
dim(residuals)

residuals_t = data.frame( t(residuals) )
dim(residuals_t)
rownames(residuals_t) = rownames(rna_seq)
colnames(residuals_t)  = colnames(rna_seq)

residuals_t_scale = f_row_scale(residuals_t)

gene = rownames(residuals_t)

gene_matrix = cbind(gene, residuals_t_scale)

write.table(gene_matrix, file = './data/462samples_peer/rnaseq/GEUVADIS.Gene.DATA_MATRIX', quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)















