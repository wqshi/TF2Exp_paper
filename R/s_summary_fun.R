source('~/R/s_function.R', chdir = T)
source('s_gene_regression_fun.R')

f_merge_data <- function(data_dir, file_pattern, subset_genes = NULL, quite = FALSE){

    feature_files = list.files(data_dir, pattern = file_pattern)
    
    cat('Total files', length(feature_files), '\n')
    feature_merge = data.frame()
    if (!is.null(subset_genes)){
        feature_files = intersect(feature_files, paste0(subset_genes, '.enet.features'))
    }
    feature_list=list()
    library(plyr)
    for (feature_file in feature_files){
        #cat(feature_file,'\n')
        result <- try(
            gene_feature<-read.table(f_p('%s/%s', data_dir, feature_file), header = TRUE),
            silent = TRUE
            
            )
        if (class(result) == "try-error"){
            if (quite == TRUE){
                print(result)
            }
            
            cat('Error in ', feature_file, '\n')
            next
        }
        feature_list[[feature_file]] = gene_feature
        #feature_merge = rbind(feature_merge, gene_feature)
    }
    feature_merge = rbind.fill(feature_list)
    return (feature_merge)
}

f_preprocess_feature_data <- function(data_dir, file_pattern, subset_genes = NULL, quite = FALSE, debug = FALSE){
    #Merge non-zeor coefficient features together.
    #Count all the features together.
    #If the merged file exists and newer than the target files, load the merged data.
    if(f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
    }
    
    feature_files = list.files(data_dir, pattern = file_pattern)
    
    cat('Total files', length(feature_files), '\n')
    feature_merge = data.frame()
    table_merge = data.frame()
    if (length(feature_files) == 0){
        flog.error('Empty feature files')
    }

    
    merge_file = f_p('%s/merge_features', data_dir)

    #check whether should redo the file reading.
    read_flag = TRUE
    if (file.exists(merge_file)){
        for (random_file in sample(x = feature_files, size = min(10, length(feature_files)))){
            if (file.info(f_p('%s/%s', data_dir, random_file))$mtime > file.info(merge_file)$mtime){
                read_flag = FALSE
                break
            }
        }

        if (read_flag == TRUE){
            flog.info('Read existing merged files')
            feature_merge =  read.table(file = (f_p('%s/merge_features', data_dir)), header = T)
              table_total_count = read.table(file = (f_p('%s/merge_feature_count', data_dir)), header = T)
              return (list(feature_merge = feature_merge, table_total_count = table_total_count))
        }
    }

    flog.info('Merge feature files')
    feature_list = list()
    table_list = list()
    for (feature_file in feature_files){
        cat(feature_file,'\n')
        result <- try(
            gene_feature<-read.table(f_p('%s/%s', data_dir, feature_file), header = TRUE),
            silent = TRUE
            
            )
        if (class(result) == "try-error"){
            if (quite == TRUE){
                print(result)
            }
            
            cat('Error in ', feature_file, '\n')
            next
        }
        
        gene_feature$rename = str_replace(gene_feature$name, '.*[|]', '')
        gene_feature$rename = str_replace(gene_feature$rename, '[.][0-9]*', '')
        gene_feature$rename = str_replace(gene_feature$rename, pattern = '[.][0-9]+$', '')
        gene_feature$rename = str_replace(gene_feature$rename, 'SNP.*', 'SNP')
        table_count = table(gene_feature$rename)
        feature_list[[feature_file]] = subset(gene_feature, score != 0)
        table_list[[feature_file]] = data.frame(feature = names(table_count), count = as.vector(table_count))
                                        #feature_merge = rbind(feature_merge, subset(gene_feature, score != 0))
        #table_merge = rbind(table_merge, data.frame(feature = names(table_count), count = as.vector(table_count)))
    }
    library(plyr)
    feature_merge = rbind.fill(feature_list)
    table_merge = rbind.fill(table_list)
    
    table_total_count <- table_merge %>% group_by(feature) %>% dplyr::summarise(total_count = sum(count))
    write.table(feature_merge, file = (f_p('%s/merge_features', data_dir)), quote = FALSE, row.names = FALSE)
    write.table(table_total_count, file = (f_p('%s/merge_feature_count', data_dir)), quote = FALSE, row.names = FALSE )
    return (list(feature_merge = feature_merge, table_total_count = table_total_count))
    
}



f_add_adjust_rsq <- function(input_performance, total_samples){
    
    input_performance$adjust.rsq = f_adjust_r_square(input_performance$performance,
                     p=input_performance$selected_features -1,
                     n = total_samples)

    input_performance$rsq = input_performance$performance
    input_performance$performance = input_performance$adjust.rsq

    input_performance$performance[input_performance$performance < 0] =0
    input_performance$performance[ input_performance$selected_features > total_samples ] =0
    return (input_performance)
}

f_summary_regression_results <- function(batch_name, chr_str, mode_name, rsync_flag = TRUE, return_features = FALSE){
    
    output_dir = f_p("./data/%s/rnaseq/%s/%s/", batch_name, chr_str, mode_name)
      
    #Rsync the features files back to the clustdell
    if (rsync_flag == TRUE){
        rsync_cmd = f_p("rsync -rav  --include '*enet*' --exclude '*' shi@clustdell.cmmt.ubc.ca:/home/shi/expression_var/data/%s/rnaseq/%s/%s/   ~/expression_var/data/%s/rnaseq/%s/%s/", batch_name, chr_str, mode_name, batch_name, chr_str, mode_name)
        system(rsync_cmd)
    }else{
        flog.info('Skip rsync...')
    }
    
    flog.info('Merge Performance')
    performance_df = f_merge_data(output_dir, '.*enet$')
    #subset_genes = as.character(subset(performance_df, performance > 0.01)$gene)
    subset_genes = performance_df$gene
    
    flog.info('Merge features')
    #features_df=f_merge_data(output_dir, '.*enet.features$', subset_genes, quite = FALSE)
    features_list = f_preprocess_feature_data(output_dir,  '.*enet.features.gz$')
    features_df = features_list$feature_merge
    
    head(performance_df)
        
    head(features_df)
    selected_feature_freq=as.data.frame(table(str_replace(features_df$name, '[|].*', '')))
    rownames(selected_feature_freq ) = selected_feature_freq$Var1
    performance_df$selected_features = selected_feature_freq[performance_df$gene,'Freq']
    #hist(performance$selected_features)
    cat('Mean performance:',mean(performance_df$performance), '\n')
    #hist(performance_df$performance)

    good_predictions <- performance_df %>% filter(performance > 0.25)
    
    cat(nrow(good_predictions), 'out of', nrow(performance_df), 'have good predictions(Rsquared > 0.25)','\n')
    write.table(features_df, f_p('%s/features', output_dir), quote =FALSE, sep = '\t', row.names= FALSE)

    write.table(performance_df, f_p('%s/performance', output_dir), quote =FALSE, sep = '\t', row.names= FALSE)
    rownames(performance_df) = performance_df$gene
    performance_df = f_add_adjust_rsq(performance_df, total_samples = 370)
    
    if (return_features == TRUE){
        return (list(features =features_df, performance = performance_df, feature_total_count = features_list$table_total_count ))
    }
    
    return (performance_df)
}



f_collect_performance_in_mode_list <- function(batch_name, chr_str, mode_list, rsync_flag){
    collected_performance = NULL
    for (loc_mode in names(mode_list)){
        cat('======', loc_mode, '========\n')
        
        result <- try(
            performance <- f_summary_regression_results( batch_name, chr_str, mode_list[loc_mode], rsync_flag),
            silent = TRUE
        )

        if (class(result) == "try-error"){
            cat('Error in ', chr, mode_list[loc_mode], '\n')
            next
        }
        
        
        #= f_adjust_r_square(performance$performance, performance$)
        performance$mode = loc_mode
        if(is.null(collected_performance)){
            collected_performance = performance
        }else{
            shared_cols = intersect(colnames(collected_performance), colnames(performance))
            collected_performance = rbind(collected_performance[, shared_cols], performance[, shared_cols])
        }
    #head10(performance)
    
    }
    return (collected_performance)
}


f_plot_two_columns <- function(input_data, col_names){
    colnames(input_data) = col_names
    input_data = as.data.frame(input_data)
    library(ggplot2)
    p<-ggplot(input_data, aes_string(x =  col_names[1], y = col_names[2])) + geom_point() + geom_abline()
    print(p)
    return (input_data)
}


f_get_gene_feature_table <- function(loc_batch_name, chr_str, mode_list, target_gene){
    
    feature_file = f_p('./data/%s/rnaseq/%s/%s/%s.enet.features', loc_batch_name, chr_str, mode_list[1], target_gene)
    feature_table = read.table(feature_file, header = T)
    feature_table$abs_score = abs(feature_table$score)
    sort_table = f_sort_by_col(feature_table, index = 'abs_score', T)
    sort_table$abs_score = NULL
    return (sort_table)
}



f_collect_performance_in_multiple_chrs<- function(loc_batch_name, mode_list, chr_num_list){    
    sample_performance_merge = data.frame()
    for (chr_num in chr_num_list){
        
        chr_str = paste0('chr', chr_num)
        dir.create(f_p('./data/%s/rnaseq/%s/', loc_batch_name, chr_str))
        
        sample_performance = f_collect_performance_in_mode_list(loc_batch_name, chr_str, mode_list, rsync_flag = TRUE)
        sample_performance = as.data.frame(sample_performance)
        
        colnames(sample_performance)
        head(sample_performance)
        sample_performance$chr = chr_str
        sample_performance_merge = rbind(sample_performance_merge, sample_performance)
        
    }

    return (sample_performance_merge)
}

f_summary_selected_features <- function(loc_batch_name, chr_list, loc_mode, output_feature_file, performance_threshold, return_tf =FALSE){
    
    accurate_features = data.frame()
    for (chr_str in chr_list){
        
        return_list = f_summary_regression_results(loc_batch_name, chr_str, loc_mode, rsync_flag = FALSE, return_features = TRUE)
        
        performance_df = return_list$performance
        features_df = return_list$features
        accurate_genes = rownames(subset(performance_df, performance > performance_threshold))

        library(plyr)
        features_df[, c('gene', 'feature')] = ldply(str_split(features_df$name,'[|]'))

        accurate_features = rbind( accurate_features, features_df[features_df$gene %in% accurate_genes,])
    }
    dim(accurate_features)

        
    tf_features=str_replace_all(accurate_features$feature, '[.][0-9]*$|[.]rs[0-9]*|promoter|enhancer|P_E|TF.overlap|esv.*', '')
    tf_features = str_replace_all(tf_features, '^[.]|[.]CEU|[.]TSI|[.]GBR|[.]YRI|[.]FIN', '')
    
    tf_features_count=sort(table(tf_features), decreasing = TRUE)
    
    feature_rank_table = sort(table(str_replace(accurate_features$feature, '[.][0-9]*$|[.]rs[0-9]*', '')), decreasing = TRUE)
    
    if(return_tf == TRUE){
        write.table(tf_features_count, output_feature_file, quote = FALSE, sep = '\t')
        return (tf_features_count)
    }else{
        write.table(feature_rank_table, output_feature_file, quote = FALSE, sep = '\t')
        return (feature_rank_table)
        
    }

    
}


f_plot_performance_and_stats_test <- function(input_performance, mode1, mode2, save_flag = FALSE ,...){

    intersect_genes = intersect(subset(input_performance, mode == mode1)$gene,
                                subset(input_performance, mode == mode2)$gene)    
    loc_performance <- input_performance%>% filter(gene %in% intersect_genes)
    mode1_perf = subset(loc_performance, mode == mode1)$performance
    mode2_perf = subset(loc_performance, mode == mode2)$performance
    if(length(mode1_perf) <3){
        flog.error('Few observations %s %s %s', length(mode1_perf), mode1, mode2)
        return (0)
    }
    
    wilcox_stats=wilcox.test(mode1_perf, mode2_perf, paired = T, alternative = 'greater')
    p <- qplot(mode1_perf, mode2_perf, alpha = I(1/3), size = I(1)) + geom_abline() + xlab(f_p("R-square of %s", mode1)) + ylab(f_p("R-square of %s", mode2))
    bigger_count=sum( mode1_perf > mode2_perf)
    less_count = sum(mode1_perf < mode2_perf)
    equal_count = sum(mode1_perf == mode2_perf)
    
    pp<-p + theme_Publication(base_size = 14) +   theme(axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3)) #+ annotate("text", label = f_p('p-value: %.1e \n %s vs %s \n %s out of %s', wilcox_stats$p.value, bigger_count, less_count, equal_count, length(mode1_perf) ),
              #       x = 0.2, y = 0.5, size = 6, colour = "red")
    
    print(pp)
}


f_venn_plot_overlap <- function(sample_performance_merge, intersect_genes, threshold){
    suppressMessages(library(limma))
    if ('TF' %in% sample_performance_merge$mode){
        TF <- subset( sample_performance_merge, gene  %in% intersect_genes & mode == 'TF')$performance >= threshold
        SNP<- subset( sample_performance_merge, gene %in% intersect_genes & mode == 'SNP')$performance >= threshold
        All<-subset( sample_performance_merge, gene %in% intersect_genes & mode == 'All')$performance >= threshold
        c3 <- cbind(TF, SNP, All)
        a <- vennCounts(c3)
        vennDiagram(a)
    }
}
f_get_feature_file_for_one_gene <- function(target_gene, mode_name){
    feature_file = f_p('./data/462samples_quantile_rmNA/rnaseq/chr22/%s/%s.enet.features', mode_name, target_gene)
    feature_table = read.table(feature_file, header = T)
    feature_table$abs_score = abs(feature_table$score)
    sort_table = f_sort_by_col(feature_table, index = 'abs_score', T)
    sort_table$abs_score = NULL
    head(sort_table)
    sort_table$chr ='chr22'
    write.table(sort_table[,c('name', 'score')], file = f_p('%s.sort', feature_file), quote = F, sep = '\t', row.names = F)
    
    rownames(sort_table) = sort_table$name
    sort_table$name = NULL
    return (sort_table)
}

f_test_performance_drop_with_ml_variance <- function(sample_performance_merge, mode1, mode2){
    intersect_genes = intersect( subset(sample_performance_merge, mode == mode1)$gene,
                                subset(sample_performance_merge, mode == mode2)$gene )

    compare_df = data.frame(m1 = subset(sample_performance_merge, mode == mode1 & gene %in% intersect_genes)$performance,
                            m1_train = subset(sample_performance_merge, mode == mode1 & gene %in% intersect_genes)$train_performance,
                            m2 = subset(sample_performance_merge, mode == mode2 & gene %in% intersect_genes)$performance,
                            m2_train = subset(sample_performance_merge, mode == mode2 & gene %in% intersect_genes)$train_performance
                            )



    print(with( subset(compare_df, m1 > 0.05), fisher.test( m1 < m2,  m1<0.5*m1_train)))
    print(with( subset(compare_df, m1 > 0.05), table( m1 < m2,  m1<0.5*m2_train)))
}



f_create_new_mode_list <- function(mode_list, target_patt, replace_patt){
    loc_mode_list = str_replace(mode_list, target_patt, replace_patt)
    names(loc_mode_list) = names(mode_list)
    return (loc_mode_list)
}



f_compare_improvment_for_two_groups <- function(group_A, group_B, good_performance, thres = 0.01){
    good_performance$ratio = good_performance$rsq/good_performance$train_performance
    snp_perf = subset(good_performance, mode == group_A )
    rownames(snp_perf) = snp_perf$gene
    all_perf = subset(good_performance, mode == group_B )
    rownames(all_perf) = all_perf$gene
    shared_genes = intersect(snp_perf$gene, all_perf$gene)
    snp_perf = snp_perf[shared_genes,]
    all_perf = all_perf[shared_genes,]
    snp_perf$group = 'Similar'
    snp_perf$group[snp_perf$performance - all_perf$performance > thres] = group_A
    snp_perf$group[snp_perf$performance - all_perf$performance < -1 * thres] = group_B

    head(snp_perf)
    snp_perf$group = factor(snp_perf$group, levels=c(group_A, group_B, 'Similar'))
    levels(snp_perf$group) = paste0(levels(snp_perf$group), ':' ,as.vector(table(snp_perf$group)))
    snp_perf$better = snp_perf$group

    print(ggplot(snp_perf, aes(train_performance, color = better)) + geom_density() + ggtitle(f_p('%s-%s test R-square', group_A, group_B)))
    print(ggplot(snp_perf, aes(ratio, color = better)) + geom_density() + xlab(f_p('%s model test/train R-square ratio', group_A)) +
          ggtitle(f_p(' Compare %s with %s model at diff threshold %s', group_A, group_B, thres )) + theme_Publication(base_size = 14) +
          theme(legend.position = 'right', legend.direction='vertical'))
                                        #ggplot(diff_df, aes(train, color = group)) + geom_density()
}


f_plot_fitting_ratio_for_two_groups <- function(group_A, group_B, good_performance){
    good_performance$ratio = good_performance$rsq/good_performance$train_performance
   
    plot_data = subset(good_performance, mode %in% c(group_A, group_B))
    shared_genes = table(as.vector(plot_data$gene) ) == 2
    shared_genes = names(shared_genes)[shared_genes]
    plot_data = subset(plot_data, gene %in% shared_genes)
    
    table(plot_data$mode)
    
    print(ggplot(plot_data, aes(ratio, color = mode)) + geom_density() + xlab(f_p('Model test/train R-square ratio')) +
          ggtitle(f_p(' Compare model fitness %s with %s model', group_A, group_B)) + theme_Publication(base_size = 14) +
          theme(legend.position = 'right', legend.direction='vertical'))
                                        #ggplot(diff_df, aes(train, color = group)) + geom_density()
}
