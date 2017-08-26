
chr_str = 'chr22'
setwd('~/expression_var/R/')
source('~/R/s_function.R', chdir = TRUE)
library(mclust)
library(mixtools)

f_read_expression_data <- function(chr_str, data_dir){
    loc_file=f_p('./data/%s/rnaseq/transcript_data.bed', data_dir, chr_str)
    loc_data = read.table(loc_file, header = TRUE, na.strings = '.')
    rownames(loc_data) = loc_data$transcript_id
    return (loc_data)
}



sailfish_data = f_read_expression_data('chr22', '462samples_sailfish')
quantile_data = f_read_expression_data('chr22', '462samples_quantile_rmNA')
peer_data = f_read_expression_data('chr22', '462samples_snyder_norm')

cat('Size of expression data:', dim(peer_data), '\n')

colnames(peer_data)

set.seed(11)

top20=read.table(f_p('./data/462samples_snyder_norm/rnaseq/chr22/top20'))
selected_genes = sample(sailfish_data$transcript_id, 50)
selected_genes =  str_replace(row.names(top20), '[.][0-9]+', '')

library(reshape2)
sample_cols = grep('(HG|NA)[0-9]', colnames(peer_data), value = TRUE)
peer_subset = melt(peer_data[selected_genes, c(sample_cols, 'transcript_id')])
peer_subset$method = 'peer'
dim(peer_data)
head(peer_subset)


sailfish_subset = melt(sailfish_data[selected_genes, c(sample_cols, 'transcript_id')])
sailfish_subset$method = 'sailfish'
head(sailfish_subset)
dim(sailfish_subset)


quantile_subset = melt(quantile_data[selected_genes, c(sample_cols, 'transcript_id')])
quantile_subset$method = 'quantile'
dim(quantile_subset)

merge_data = rbind(sailfish_subset, quantile_subset, peer_subset)

head(merge_data)
loc_gene = 'ENSG00000187642'

library(dplyr)

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

for (loc_gene in unique(merge_data$transcript_id)){

    
    loc_plot <- ggplot(subset(merge_data, transcript_id == loc_gene), aes(x=value)) + geom_density() + facet_wrap(~method, scales = "free")
    ggsave(f_p('./data/r_results/gene_distrbution/%s.tiff', loc_gene), plot = loc_plot, dpi = 75)

    gene_data = subset(merge_data, transcript_id == loc_gene & method == 'peer')
    dens = densityMclust(gene_data$value, modelName='V')
    #summary(dens)
    #summary(dens, parameters = TRUE)
    #plot(dens, what = "BIC")

    
    #plotDensityMclust1(dens, data=gene_data$value, hist.border = "grey")
    
    

   gene_data = f_sort_by_col(gene_data, index = 'value')
   hist_plot<-ggplot(gene_data) + geom_histogram(aes(value, ..density..), binwidth = 0.2, colour = "black", fill = "white")
  
    
    for (i in 1:dens$parameters$variance$G ){
        hist_plot <- hist_plot +   stat_function(geom = "line", fun = plot_mix_comps,
                arg = list(dens$parameters$mean[i], dens$parameters$variance$sigmasq[i], dens$parameters$pro[i]),
                colour = "red", lwd = 1)
    }
    #i = 1
    #plot_mix_comps(sort(gene_data$value), dens$parameters$mean[i], dens$parameters$variance$sigmasq[i], dens$parameters$pro[i]  )
    
    
    density_plot<- hist_plot + ggtitle(f_p('\n\n%s: %s', loc_gene, dens$G))
    #density_plot
    
    ggsave(f_p('./data/r_results/gene_distrbution/%s.cluster.tiff', loc_gene), plot = density_plot, dpi = 75)
    
    
    
}













