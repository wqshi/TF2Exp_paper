library(peer)
source('~/R/s_function.R', chdir = TRUE)
#install.packages('preprocessCore')
#source("http://bioconductor.org/biocLite.R")
#biocLite("preprocessCore")
library(preprocessCore)


f_row_scale <- function(df){
    return (t(scale(t(df))))
}


batch_name = '462samples_sailfish'
#batch_name = '54samples_evalue'
batch_output_dir = f_p('./data/%s/', batch_name)

rna_seq_raw = read.csv(file = f_p('%s/rnaseq/GEUVADIS.Gene.DATA_MATRIX', batch_output_dir) , sep=' ', header =TRUE)
rownames(rna_seq_raw) = rna_seq_raw$gene

rna_seq_raw$gene = NULL
rna_seq_raw[is.na(rna_seq_raw)] = 0

head(rna_seq_raw)

library(reshape2)

set.seed(11)
#check the assumption for quantile normalization
random_samples = sample(colnames(rna_seq_raw), size = 15)
long_exp = melt(rna_seq_raw[, random_samples])
head(long_exp)
sample_plot <- ggplot(long_exp, aes(x = value, group = variable, color = variable)) + geom_density() + xlim(c(0,50))

#Check the distribution for genes.
random_genes = sample(rownames(rna_seq_raw), size = 20)
long_gene_exp = melt(cbind(gene = random_genes, rna_seq_raw[random_genes,]))
head(long_gene_exp)
#This is useless
gene_plot <- ggplot(long_gene_exp, aes(x = value, group = gene, color = gene)) + geom_density() + xlim(c(0,10))





##Quantile normalization

rna_seq_quantile = normalize.quantiles(as.matrix(rna_seq_raw))
colnames(rna_seq_quantile) = colnames(rna_seq_raw)

head(rna_seq_quantile[,1:10])
head(rna_seq_raw[,1:10])

rna_seq = scale(rna_seq_raw)

rna_seq_scale = f_row_scale(rna_seq_quantile)

head(rna_seq_scale[,1:10])

dim(rna_seq_scale)


gene = data.frame( gene = rownames(rna_seq_raw) )
write_data = cbind(gene, rna_seq_scale)
#head(write_data[,1:10])


write.table(write_data, file = './data/462samples_quantile/rnaseq/GEUVADIS.Gene.DATA_MATRIX', quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)

##########variance stablize and quantile normalization
peaksMat = as.matrix(rna_seq_raw)
peaksMat_asinh = asinh(peaksMat) 

asinh_quantile = normalize.quantiles(as.matrix(peaksMat_asinh))
colnames(asinh_quantile) = colnames(rna_seq_raw)


asinh_scale = f_row_scale(asinh_quantile)



gene = data.frame( gene = rownames(rna_seq_raw) )
asinh_write_data = cbind(gene, asinh_scale)
#head(write_data[,1:10])

write.table(asinh_write_data, file = './data/462samples_var_quantile/rnaseq/GEUVADIS.Gene.DATA_MATRIX', quote = FALSE, sep = ' ', row.names = FALSE, col.names = TRUE)






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















