
loc_batch_name = '358samples_regionkeepLow'
collection_name = 'peer358corRmdup'    



results_dir = f_p('./data/%s/rnaseq/results/%s/', loc_batch_name, collection_name)
figure_dir = f_p('./data/%s/rnaseq/results/%s/figures', loc_batch_name, collection_name)
dir.create(results_dir)

fdr_threshold = 0.05

CTCF_data = read.table('./data/raw_data/CEU47/CTCF_validation_results', header = T)
PU1_data = read.table('./data/raw_data/CEU47/PU1_validation_results', header = T)

combined_data = rbind(CTCF_data, PU1_data)


combined_data$sig_tf_gene = combined_data$tf_gene_fdr < fdr_threshold
combined_data$sig_tf_deepsea = combined_data$tf_deepsea_fdr < fdr_threshold

head(combined_data$sig_tf_gene)

ggplot(combined_data, aes(tf_gene_coef, gene_deepsea_coef, colour = sig_tf_gene, size = sig_tf_deepsea ))+
    geom_point(alpha = 0.5) + geom_abline() + coord_fixed() + xlim(c(-0.8,0.8)) + ylim(c(-0.8,0.8)) +
    xlab('Correlation between real TF binding and gene expression') + ylab('Correlation between DeepSEA and gene expression') +
    scale_color_discrete( labels = c("Significant TF-gene correlation", "No significance"), guide = guide_legend(title = 'Size')) +
    scale_size_discrete( labels = c("Significant TF-DeepSEA correlation", "No significance"), guide = guide_legend(title = 'Color')) +
    theme_Publication(12)
    #scale_size_discrete( labels = c("No significance", "Significant TF-Expression correlation"), guide = guide_legend(title = NULL))
    

combined_data$sig_tf_gene = 'Sig TF-Gene'
combined_data$sig_tf_deepsea = 'Sig TF-DeepSEA'

combined_data$sig_gene_deepsea = 'Sig DeepSEA-Gene'

combined_data$sig_tf_gene[combined_data$tf_gene_fdr > fdr_threshold] = 'No Sig'
combined_data$sig_tf_deepsea[combined_data$tf_deepsea_fdr > fdr_threshold] = 'No Sig'

combined_data$sig_gene_deepsea[combined_data$gene_deepsea_fdr > fdr_threshold] = 'No Sig'

combined_data$Color = factor(combined_data$sig_tf_gene, levels = c('Sig TF-Gene', 'No Sig'))
#combined_data$Fill = factor(combined_data$sig_tf_deepsea, levels = c('Sig TF-DeepSEA', 'No Sig'))
combined_data$Fill = factor(combined_data$sig_gene_deepsea, levels = c('Sig DeepSEA-Gene', 'No Sig'))

table(combined_data$Fill)

sig_data = subset(combined_data, tf_gene_fdr < fdr_threshold)
cor_obj = cor.test( sig_data$gene_deepsea_coef, sig_data$tf_gene_coef, paired = T)

cor_obj
table(combined_data$factor)
table(combined_data$tf_gene_fdr < fdr_threshold)

table(combined_data$Fill)

ggplot(combined_data, aes(tf_gene_coef, gene_deepsea_coef, colour = Color)) +
    theme_Publication(12) + geom_point(size = 2) + theme(legend.position = 'right', legend.direction = 'vertical') +
    scale_colour_manual(values = c("#F8766D", "grey"))  + geom_abline() + coord_fixed() + xlim(c(-0.8,0.8)) + ylim(c(-0.8,0.8)) +
    xlab('Correlation coeficient b/w TF and gene expression') +
    ylab('Correlation coeficient b/w DeepSEA and gene expression') +
    annotate("text", x = -0.6, y = 0.7, label = f_p("P: %.2e", cor_obj$p.value)) +
    annotate("text", x = -0.6, y = 0.62, label = f_p("Corr. : %.2f", cor_obj$estimate))

ggsave(f_p('%s/validation_figure.tiff', figure_dir), width = 7, height =7)

write.table(table(combined_data$tf_gene_fdr < fdr_threshold, combined_data$factor), file = f_p('%s/validation_results.txt', results_dir), quote = F)
write('##TRUE means sig TF-gene correlation', file = f_p('%s/validation_results.txt', results_dir),  append = T)
 
ggplot(combined_data, aes(tf_deepsea_coef, tf_deepsea_fdr)) + geom_point()


range(abs(subset(combined_data, tf_deepsea_fdr < 0.05)$tf_deepsea_coef))


