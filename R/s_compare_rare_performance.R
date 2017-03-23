source('~/R/s_ggplot2_theme.R', chdir = TRUE)
source('s_summary_fun.R')
#source('./s_project_data.R')
source('s_project_funcs.R')

loc_batch_name = '358samples_regionkeepLow'
collection_name = 'peer358corRmdup'
figure_dir = f_p('./data/%s/rnaseq/results/%s/figures', loc_batch_name, collection_name)
dir.create(figure_dir)
results_dir = f_p('./data/%s/rnaseq/results/%s/', loc_batch_name, collection_name)

rare_number_file = f_p('%s/paper_rare_numbers.txt', results_dir)




mode_list = modes_list[[collection_name]][c('TF', 'SNP','All', 'SNPinTF', 'random' ,'TFaddInteract','TFsnpMatch', 'TFaddPenalty')]
chr_num_list = c(1:22)

all_performance_merge = f_collect_performance_in_multiple_chrs('358samples_regionkeepLow', mode_list[c('TF')], chr_num_list, debug = F)
snpOnly_performance_merge = f_collect_performance_in_multiple_chrs('358samples_regionkeepLow_snpOnly', mode_list['TF'], chr_num_list)
rare_performance_merge = f_collect_performance_in_multiple_chrs('358samples_regionkeepLow_rareVar', mode_list[c('TF')], chr_num_list, debug = F)

all_performance_merge$mode_full = paste0(all_performance_merge$mode, '.All')
snpOnly_performance_merge$mode_full = paste0(snpOnly_performance_merge$mode, '.SNP')
rare_performance_merge$mode_full = paste0(rare_performance_merge$mode, '.rare')

colnames(all_performance_merge)
colnames(snpOnly_performance_merge)



dim(all_performance_merge)

shared_cols = intersect( intersect(colnames(rare_performance_merge), colnames(all_performance_merge)), colnames(snpOnly_performance_merge) )



performance_merge = rbind(all_performance_merge[, shared_cols], snpOnly_performance_merge[,shared_cols], rare_performance_merge[, shared_cols])
table(performance_merge$mode_full, performance_merge$chr)
with(all_performance_merge, table(mode_full, chr))

performance_sort <- performance_merge %>% arrange(mode_full, performance) %>%
    group_by(mode_full) %>% dplyr::mutate( rank = row_number(mode_full)/length(mode_full) )



random_tf_plot <-
    ggplot(performance_sort, aes(rank, performance, color = mode_full)) + geom_point(size = 0.1) +
    theme_Publication( base_size = 12) + 
    theme(legend.position = c(0.2,0.8), legend.direction = 'vertical') + geom_hline(yintercept=0.05, linetype = '1F') +
    xlab('Normalized rank of investigated genes') + ylab('Model Performance (R2)') +
    scale_y_continuous(breaks=c(0,0.05,0.2,0.4,0.6,0.8)) +
    scale_color_discrete( guide = guide_legend(override.aes = list(size=2), title = NULL))

random_tf_plot

ggsave(filename = f_p('%s/perf_rare_snp_all.tiff', figure_dir), plot = random_tf_plot, width = 7, height = 7, units = 'in')

library(tidyr)
head(performance_merge)



##Output: mean performance of the rare.
table(performance_merge$mode_full)
performance_merge %>% group_by(mode_full) %>% dplyr::summarise(mean_perf = mean(performance, na.rm = T))








plot_data <- performance_merge %>% dplyr::select(gene, performance, mode_full) %>%
    spread(mode_full, performance) %>% mutate(perf_diff = TF.All - TF.SNP)

##Output: how much is improved
plot_data %>% dplyr::select(gene, perf_diff) %>% dplyr::summarise(sum(perf_diff>0, na.rm = T)/sum(!is.na(perf_diff)))


cor_test = cor.test(plot_data$TF.rare, plot_data$perf_diff)
cor_test$p.value = f_get_exact_pvalue(cor_test)

f_write_paper_results('======s_compare_rare_performance==========', date(), rare_number_file)
cat(collection_name, loc_batch_name, '\n', file = rare_number_file, append = TRUE)
f_write_paper_results('Test the correlation between rare improvment and rare performance', '', rare_number_file)
f_write_paper_results('Pearson correlation p-value', cor_test$p.value, rare_number_file)
f_write_paper_results('coefficient', cor_test$estimate, rare_number_file)
#TFsnpMatch_test = cor.test(plot_data$TFsnpMatch.rare, plot_data$perf_diff)
#TFsnpMatch_test$p.value = f_get_exact_pvalue(TFsnpMatch_test)



#Output figure:
rare_improvement <- ggplot(plot_data, aes(TF.rare, perf_diff)) + geom_point(alpha = 0.4) +
    xlab('R2 of uncommon variants models') +
    ylab('R2 improvement by adding uncommon variants on top of common variants') + 
    annotate(geom = 'text', x = 0.05, y = 0.5, label = f_p('P-value: %.1e\nCorr: %.2f     ', cor_test$p.value, cor_test$estimate), hjust = 0) +
    theme_Publication(12) + #geom_smooth() +
    geom_hline(yintercept=0, color = 'grey' ) 

rare_improvement


improvement_without_line <- ggplot(plot_data, aes(TF.rare, perf_diff)) + geom_point(alpha = 0.4) +
    xlab('Performance of uncommon variants models') +
    ylab('Performance improvement by adding uncommon variants \n on top of common variants') + 
    #annotate(geom = 'text', x = 0.05, y = 0.5, label = f_p('P-value: %.1e\nCorr: %.2f     ', cor_test$p.value, cor_test$estimate), hjust = 0) +
    theme_Publication(12) + #geom_smooth() +
    geom_hline(yintercept=0, color = 'grey' )

improvement_without_line

ggsave(filename = f_p('%s/rare_improvement.tiff', figure_dir), plot = improvement_without_line, width = 7, height = 7, units = 'in')

