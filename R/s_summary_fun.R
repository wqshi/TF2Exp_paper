
source('~/R/s_function.R', chdir = T)
source('s_gene_regression_fun.R')
source('s_project_funcs.R')
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


f_rename_features <- function(feature_merge){

    feature_merge$rename = str_replace(feature_merge$name, '.*[|]', '')

    feature_merge$rename = gsub(pattern = '(?<!PU|USF|EGR)[.][0-9]+$', '', feature_merge$rename, perl = T)

    
    feature_merge$rename = gsub(pattern = '(?<=promoter|enhancer)[.]', '', feature_merge$rename, perl = T)
    
    feature_merge$rename = str_replace(feature_merge$rename, 'SNP.*', 'SNP')
    
    return (feature_merge$rename)
}



f_preprocess_feature_data <- function(data_dir, file_pattern, subset_genes = NULL, quite = FALSE, read_flag = TRUE, debug = FALSE){
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

    merge_file = f_p('%s/merge_features.control', data_dir)

    #check whether should redo the file reading.
    
    if (file.exists(merge_file) & read_flag == TRUE){
        for (random_file in sample(x = feature_files, size = min(10, length(feature_files)))){
            if (file.info(f_p('%s/%s', data_dir, random_file))$mtime > file.info(merge_file)$mtime){
                read_flag = FALSE
                break
            }
        }
        
        if (read_flag == TRUE){
            flog.info('Read existing merged files')
            feature_merge =  read.table(file = (f_p('%s/merge_features', data_dir)), header = T)

            if (file.info(f_p('%s/merge_features.control', data_dir))$size > 10){
                control_merge =  read.table(file = (f_p('%s/merge_features.control', data_dir)), header = T)
                control_merge$rename = f_rename_features(control_merge)
            }else{
                control_merge = NULL
            }
            table_total_count = read.table(file = (f_p('%s/merge_feature_count', data_dir)), header = T)

            feature_merge$rename = f_rename_features(feature_merge)
            
            return (list(feature_merge = feature_merge, table_total_count = table_total_count, control_merge = control_merge))
        }
    }#Else: do the following
    feature_file = feature_files[1]
    flog.info('Merge feature files')
    feature_list = list()
    control_list = list()
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
        #gene_feature$rename = str_replace(gene_feature$rename, '[.][0-9]*', '')
        #str_replace(gene_feature$rename, pattern = '[.][0-9]+$', '')
        gene_feature$rename = gsub(pattern = '(?<!PU|USF|EGR)[.][0-9]+$', '', gene_feature$rename, perl = T)

        gene_feature$rename = str_replace(gene_feature$rename, 'SNP.*', 'SNP')
        table_count = table(gene_feature$rename)
        feature_list[[feature_file]] = subset(gene_feature, score != 0)
        control_features=gene_feature[!grepl('Intercept', gene_feature$name) & grepl('promoter|enhancer', gene_feature$name) & gene_feature$score == 0,]
        target_size =                    sum(!grepl('Intercept', gene_feature$name) & grepl('promoter|enhancer', gene_feature$name) & gene_feature$score != 0)
        control_sel=sample(x=1:nrow(control_features), size = target_size, replace = target_size > nrow(control_features) )

        if (length(control_sel) > 0){
            control_list[[feature_file]]  = control_features[control_sel,]
        }
        
        table_list[[feature_file]] = data.frame(feature = names(table_count), count = as.vector(table_count))
                                        #feature_merge = rbind(feature_merge, subset(gene_feature, score != 0))
        #table_merge = rbind(table_merge, data.frame(feature = names(table_count), count = as.vector(table_count)))
    }
    library(plyr)
    
    feature_merge = rbind.fill(feature_list)
    control_merge = rbind.fill(control_list)
    table_merge = rbind.fill(table_list)
    library(dplyr)
    table_total_count <- table_merge %>% dplyr::group_by(feature) %>% dplyr::summarise(total_count = sum(count))
    write.table(feature_merge, file = (f_p('%s/merge_features', data_dir)), quote = FALSE, row.names = FALSE)
    write.table(control_merge, file = (f_p('%s/merge_features.control', data_dir)), quote = FALSE, row.names = FALSE)
    write.table(table_total_count, file = (f_p('%s/merge_feature_count', data_dir)), quote = FALSE, row.names = FALSE )
    return (list(feature_merge = feature_merge, table_total_count = table_total_count, control_merge = control_merge))
    
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

f_summary_regression_results <- function(batch_name, chr_str, mode_name, rsync_flag = TRUE, return_features = FALSE, read_flag = TRUE, debug = FALSE){

    if(f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        
    }
    
    #read_flag: read the compiled data.
    output_dir = f_p("./data/%s/rnaseq/%s/%s/", batch_name, chr_str, mode_name)
        
    #Rsync the features files back to the clustdell
    if (rsync_flag == TRUE){        
        delete_cmd = f_p('find ~/expression_var/data/%s/rnaseq/%s/%s/ -name *enet* -delete', batch_name, chr_str, mode_name)
        #rsync_cmd = f_p("rsync -rav  --include '*enet*' --exclude '*' shi@clustdell.cmmt.ubc.ca:/home/shi/expression_var/data/%s/rnaseq/%s/%s/   ~/expression_var/data/%s/rnaseq/%s/%s/", batch_name, chr_str, mode_name, batch_name, chr_str, mode_name)
        rsync_cmd = f_p("rsync -rav  --include '*' shi@clustdell.cmmt.ubc.ca:/home/shi/expression_var/data/%s/rnaseq/%s/%s/   ~/expression_var/data/%s/rnaseq/%s/%s/", batch_name, chr_str, mode_name, batch_name, chr_str, mode_name)
        #system(delete_cmd)
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
    features_list = f_preprocess_feature_data(output_dir,  '.*enet.features.gz$', read_flag = read_flag, debug = F)
    features_df = features_list$feature_merge
    features_df$gene = str_replace(features_df$name, '[|].*', '')

    #f_ASSERT( setequal(performance_df$gene, features_df$gene), 'Unequal sets of genes in features and performance' )
    
    
    head(performance_df)
        
    head(features_df)
    selected_feature_freq=as.data.frame(table(str_replace(features_df$name, '[|].*', '')))
    rownames(selected_feature_freq ) = selected_feature_freq$Var1
    performance_df$selected_features = selected_feature_freq[as.character(performance_df$gene),'Freq']
    #hist(performance$selected_features)
    cat('Mean performance:',mean(performance_df$performance), '\n')
    #hist(performance_df$performance)
    library(dplyr)
    good_predictions <- performance_df %>% filter(performance > 0.25)
    
    cat(nrow(good_predictions), 'out of', nrow(performance_df), 'have good predictions(Rsquared > 0.25)','\n')
    write.table(features_df, f_p('%s/features', output_dir), quote =FALSE, sep = '\t', row.names= FALSE)

    write.table(performance_df, f_p('%s/performance', output_dir), quote =FALSE, sep = '\t', row.names= FALSE)
    rownames(performance_df) = performance_df$gene
    performance_df_adjust = f_add_adjust_rsq(performance_df, total_samples = 370)
    
    if (return_features == TRUE){
        return (list(features =features_df, performance = performance_df_adjust, feature_total_count = features_list$table_total_count, control = features_list$control_merge ))
    }
    
    return (performance_df)
}



f_collect_performance_in_mode_list <- function(batch_name, chr_str, mode_list, rsync_flag, debug = F){

    if(f_judge_debug(debug)){

         ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
        
    }
    
    collected_performance = NULL
    for (loc_mode in names(mode_list)){
        cat('======', loc_mode, '========\n')

        
        result <- try(
            performance <- f_summary_regression_results( batch_name, chr_str, mode_list[loc_mode], rsync_flag, debug = F),
            silent = TRUE
        )

        if (class(result) == "try-error"){
            cat('Error in ', chr_str, mode_list[loc_mode], '\n')
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



f_collect_performance_in_multiple_chrs<- function(loc_batch_name, mode_list, chr_num_list, debug = F){    
    if(f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        
    }

    sample_performance_merge = data.frame()
    for (chr_num in chr_num_list){
        
        chr_str = paste0('chr', chr_num)
        dir.create(f_p('./data/%s/rnaseq/%s/', loc_batch_name, chr_str))
                
        sample_performance = f_collect_performance_in_mode_list(loc_batch_name, chr_str, mode_list, rsync_flag = TRUE)
        sample_performance = as.data.frame(sample_performance)
        
        colnames(sample_performance)
        head(sample_performance)
        sample_performance$chr = chr_str
        if (nrow(sample_performance_merge) == 0){
            shared_cols = colnames(sample_performance)
            shared_cols = grep('X[0-9]', shared_cols, invert = T, value =T)
            sample_performance_merge = rbind(sample_performance_merge, sample_performance[, shared_cols])
        }else{
            shared_cols = intersect(colnames(sample_performance_merge), colnames(sample_performance))
            sample_performance_merge = rbind(sample_performance_merge[, shared_cols], sample_performance[, shared_cols])
        }
        
        
        
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


f_plot_performance_and_stats_test <- function(input_performance, mode1, mode2, publication = FALSE, save_flag = FALSE , test_stats = F, debug = F,...){


    if (f_judge_debug(debug)){

        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        
    }
    
    intersect_genes = intersect(subset(input_performance, mode == mode1)$gene,
                                subset(input_performance, mode == mode2)$gene)    
    loc_performance <- input_performance %>% filter(gene %in% intersect_genes)
    mode1_perf = subset(loc_performance, mode == mode1)$performance
    mode2_perf = subset(loc_performance, mode == mode2)$performance
    if(length(mode1_perf) <3){
        flog.error('Few observations %s %s %s', length(mode1_perf), mode1, mode2)
        return (0)
    }
    
    wilcox_stats=wilcox.test(mode1_perf, mode2_perf, paired = T, alternative = 'two.sided', conf.int =T)
    p <- qplot(mode1_perf, mode2_perf, alpha = I(1/3), size = I(1), color = I('skyblue4')) + geom_abline() + xlab(f_p("R2 of %s", mode1)) + ylab(f_p("R2 of %s", mode2))
    bigger_count=sum( mode1_perf > mode2_perf)
    less_count = sum(mode1_perf < mode2_perf)
    equal_count = sum(mode1_perf == mode2_perf)
    pp<-p + theme_Publication(base_size = 12) +   theme(axis.line.x = element_line(color="black", size = 0.3),
                                                            axis.line.y = element_line(color="black", size = 0.3))
        
    mean_diff = mean(mode1_perf) - mean(mode2_perf)
    if (publication == FALSE){
 
         ppp <- pp + annotate("text", label = f_p('P-value: %.1e \n %s vs %s \n %s out of %s', wilcox_stats$p.value, bigger_count, less_count, equal_count, length(mode1_perf) ),
                     x = 0.2, y = 0.5, size = 6, colour = "red") + annotate('text', label = f_p('X: %.2e, Y:%.2e, Mean diff %.2e', mean(mode1_perf), mean(mode2_perf), mean_diff), x = 0.3, y = 0.6) +
             annotate('text', label = f_p('X: %.2e, Y:%.2e, Median diff %.2e', median(mode1_perf), median(mode2_perf), median(mode1_perf) - median(mode2_perf)  ), x = 0.3, y = 0.65)
    }else if(publication == TRUE){
        ppp <- pp + annotate("text", label = f_p('P-value: %.1e ', wilcox_stats$p.value),
                     x = 0.2, y = 0.75, size = 5, colour = 'black')
    }else{
        ppp <- pp + annotate("text", label = f_p('P-value: %.2f ', wilcox_stats$p.value),
                     x = 0.2, y = 0.75, size = 5, colour = 'black')
        
    }

    if (test_stats == TRUE){
        ppp <- list( pvalue = wilcox_stats$p.value, median_difference = wilcox_stats$estimate,
                    mean_diff = mean(mode1_perf) - mean(mode2_perf), mean1 = mean(mode1_perf), mean2 = mean(mode2_perf))
        
    }
    
    #print(pp)

    return (ppp)
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



f_compare_improvment_for_two_groups <- function(group_A, group_B, good_performance, thres = 0.01, return_flag = FALSE, debug = FALSE){

    if(f_judge_debug(debug)){

        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        
    }
    
    good_performance$ratio = good_performance$rsq/good_performance$train_performance
    cat('Range of ratio', min(good_performance$ratio), max(good_performance$ratio), '\n')
    snp_perf = subset(good_performance, mode == group_A )
    rownames(snp_perf) = snp_perf$gene
    all_perf = subset(good_performance, mode == group_B )
    rownames(all_perf) = all_perf$gene
    shared_genes = intersect(snp_perf$gene, all_perf$gene)
    snp_perf = snp_perf[shared_genes,]
    all_perf = all_perf[shared_genes,]

    cat('Mean performance', group_A, mean(snp_perf$performance), group_B, mean(all_perf$performance), '\n')

    snp_perf$group = 'Similar'
    snp_perf$group[snp_perf$performance - all_perf$performance > thres] = group_A
    snp_perf$group[snp_perf$performance - all_perf$performance < -1 * thres] = group_B
    snp_perf$perf_diff = all_perf$performance - snp_perf$performance
        
    head(snp_perf)
    snp_perf$group = factor(snp_perf$group, levels=c(group_A, group_B, 'Similar'))
    levels(snp_perf$group) = paste0(levels(snp_perf$group), ':' ,as.vector(table(snp_perf$group)))
    snp_perf$better = snp_perf$group
    
    
    
    #print(ggplot(snp_perf, aes(train_performance, color = better)) + geom_density() + ggtitle(f_p('%s-%s test R-square', group_A, group_B)))
    density_plot <-ggplot(snp_perf, aes(ratio, color = better)) + geom_density() + xlab(f_break_long_line_into_multiple(f_p('%s model test/train R2 ratio', group_A))) +
          ggtitle(f_break_long_line_into_multiple(f_p(' Compare %s with %s at diff threshold %s', group_A, group_B, thres ))) + theme_Publication(base_size = 14) +
          theme(legend.position = 'right', legend.direction='vertical')
                                        #ggplot(diff_df, aes(train, color = group)) + geom_density()
    table(snp_perf$group)
    result <- try(
        test_obj <- wilcox.test(subset(snp_perf, grepl(group_A, better))$ratio, subset(snp_perf, grepl(group_B, better))$ratio),
        silent = TRUE
    )
    if (class(result) == "try-error"){
        print(result)
    }else{
    
        cat('Fitness difference between two groups, test p-value:', test_obj$p.value, '\n')
    }
    overall_density_plot <-ggplot(snp_perf, aes(ratio)) + geom_density() +
        xlab(f_break_long_line_into_multiple(f_p('%s model test/train R2 ratio', group_A))) +
        theme_Publication(base_size = 14)
    
    diff_plot <-ggplot(snp_perf, aes(ratio, perf_diff)) + geom_point() +
        ylab(f_break_long_line_into_multiple(f_p('%s - %s', group_B, group_A))) +
        xlab(f_break_long_line_into_multiple(f_p('%s test/train R2 ratio', group_A))) +
    ggtitle(f_break_long_line_into_multiple(f_p('Performance difference')))
    snp_perf$other_selected = all_perf$selected_features

    if (return_flag){
        return (snp_perf[,c('performance', 'SD','perf_diff', 'ratio', 'num_feature', 'selected_features', 'other_selected')])
    }else{
        return (list( diff =diff_plot, density = density_plot, overall_density = overall_density_plot ))
    }
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


modes_list = list()

modes_list$addPenalty=c(
All = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.All_other.info.normCor',
#AllfilterMinor = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.AllfilterMinor_other.info.normCor',
#AlltopTF = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.AlltopTF_other.info.normCor',
#All = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.All_other.info.normCor',
SNPinTF = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.SNPinTF_other.info.normCor',
SNP = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.SNP_other.info.normCor',
random = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.random_other.info.normCor',
TFsnpMatch = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.TFsnpMatch_other.info.normCor',
TF = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.TF_other.info.normCor',
TFaddPenalty = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.TFaddPenalty_other.info.normCor',
#,
#AllrmHic='rm.histone_model.cv.glmnet_rm.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.All_other.info.normCor',
#AllnoInteract = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.AllnoInteract_other.info.normCor',
##InterOnlySNPinTF='rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.InterOnlySNPinTF_other.info.normCor',
TFaddInteract ='rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.TFaddInteract_other.info.normCor',
TFfilterMinor ='rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.TFfilterMinor_other.info.normCor'
#AlltfShuffle = 'rm.histone_model.cv.glmnet_add.penalty_population.all_new.batch.445samples.snyder.original_batch.mode.AlltfShuffle_other.info.normCor'
)

modes_list$snpOnly=modes_list$addPenalty#[c('All', 'TF', 'SNP')]
modes_list$addPenaltyRmCor = f_create_new_mode_list(modes_list$addPenalty, 'normCor', 'rmCor')
modes_list$snyder.norm = f_create_new_mode_list(modes_list$addPenalty, 'snyder.original', 'snyder.norm')
modes_list$maxit = f_create_new_mode_list(modes_list$addPenalty, 'normCor', 'maxit1M')
modes_list$lm = f_create_new_mode_list(modes_list$addPenalty, 'normCor', 'lm')

modes_list$non_pop = f_create_new_mode_list(modes_list$snyder.norm, 'population.all', 'population.None')
modes_list$old2TFinter = f_create_new_mode_list(modes_list$non_pop, 'normCor', 'old2TFinter')
modes_list$tradR2 = f_create_new_mode_list(modes_list$non_pop, 'normCor', 'tradR2')
modes_list$tradR2keepZero = f_create_new_mode_list(modes_list$non_pop, 'normCor', 'tradR2keepZero')
modes_list$rm.penalty = f_create_new_mode_list(modes_list$tradR2keepZero, 'add.penalty', 'rm.penalty')
modes_list$rm.YRI = f_create_new_mode_list(modes_list$tradR2keepZero, 'add.penalty', 'rm.penalty_rm.YRI')
modes_list$peer = f_create_new_mode_list(modes_list$rm.YRI, 'tradR2keepZero', 'tradR2keepZeroPopPeer')
modes_list$peer = f_create_new_mode_list(modes_list$peer, 'snyder.norm', 'peer')
modes_list$peer358 = f_create_new_mode_list(modes_list$peer, '445', '358')
modes_list$peer358cor = f_create_new_mode_list(modes_list$peer358, 'trad', 'cor')

modes_list$peer358corRmdup = f_create_new_mode_list(modes_list$peer358cor, 'Peer', 'PeerRmdup')
modes_list$peer358corRmdup[['SNP']] = modes_list$peer358cor[['SNP']]
modes_list$peer358corRmdup[['SNPinTF']] = modes_list$peer358cor[['SNPinTF']]

modes_list$validation[['CTCF']] = paste0(modes_list$peer358corRmdup[['TF']], 'CTCF')
modes_list$validation[['PU1']] = paste0(modes_list$peer358corRmdup[['TF']], 'PU1')


modes_list$gtex = f_create_new_mode_list(modes_list$peer358cor, 'Peer', 'GTexRmdup')
modes_list$gtex = f_create_new_mode_list(modes_list$gtex, 'peer', 'gtex.norm')
modes_list$elastic = f_create_new_mode_list(modes_list$peer358cor, 'Rmdup', 'RmdupElastic')

modes_list$rmZero = f_create_new_mode_list(modes_list$peer358corRmdup, 'corR2keepZeroPopPeerRmdup', 'corR2PopPeerRmdup')
modes_list$rmZero[['SNP']] = modes_list$peer358corRmdup[['SNP']]
modes_list$rmZero[['All']] = modes_list$peer358corRmdup[['All']]


f_get_individaul_genotype_from_vcf_file <- function(individual_id, project_dir, batch_name){
    indiv_vcf = f_p('%s/data/%s/chr_vcf_files/chr22/%s.vcf.gz', project_dir, batch_name, individual_id)

    vcf_df = read.table(indiv_vcf)
    colnames(vcf_df) = c('chr', 'pos', 'genotype', 'ref', 'alt')

    vcf_df$add_score = 0.5
    vcf_df[vcf_df$genotype == '1|1','add_score'] = 1

    head(vcf_df)

    head(var_feature_bed_subset)

    vcf_df$snp_id = paste0(vcf_df$pos,':', vcf_df$ref, ':', vcf_df$alt)
    rownames(vcf_df) = vcf_df$snp_id
    return (vcf_df)
}


f_test_preprocess_for_one_gene <- function(target_gene, chr_str, batch_name, return_list_tf, var_feature_bed_subset, debug = FALSE){

    if (f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
   
    }
    
    ##Read the final processed data.
    test_gene <- GENE(data = data.frame(), gene_name = target_gene, chr_str = chr_str, batch_name = batch_name)
    test_gene$read_data()
    test_gene$subset_features_to_snp_contain_region('TF', debug = F)
    test_data_flag=test_gene$read_test_data(test_batch = '800samples', debug = FALSE)
    dim(test_gene$data)
    test_gene$hic_within_1Mb()
    transcript_data = test_gene$data
    rownames(transcript_data) = make.names(paste0(transcript_data$type, '-',transcript_data$feature), unique = T)
    dim(transcript_data)
    f_input_stats(t(transcript_data), batch_mode = 'TF')
    set.seed(11)
    ##Get individual genotype data
    selected_samples = sample( grep('t_', test_gene$get_samples(), value = T, invert = T), size = 10)

    individual_id = 'HG00271'
    project_dir = './'
    for (individual_id in selected_samples){
        vcf_df <-f_get_individaul_genotype_from_vcf_file(individual_id, project_dir, batch_name)
        table(vcf_df$genotype)

        ##Select one key features of the target gene
        key_features_df <- return_list_tf$features %>% separate(name, into =c('gene', 'feature'), sep ='[|]') %>% filter(gene == target_gene)
        key_features = grep('Intercept',key_features_df$feature,value = T, invert = T)
        if (length(key_features) != 1){
            key_features = key_features[1]
        }
        key_features_df[key_features_df$feature == key_features,]

        key_features_df[key_features_df$feature == key_features, c('feature_start', 'feature_end')]
        transcript_data[key_features,c('feature_start', 'feature_end')]

        feature_start = transcript_data[key_features,'feature_start']
        feature_end = transcript_data[key_features,'feature_end']
        key_tfs = transcript_data[key_features,'feature']
        key_tfs
        ##Recalculate the variants impact in the selected feature
        individual_snps_ids = (vcf_df %>% filter(pos >= feature_start, pos <= feature_end))
        all_var_in_tf <- var_feature_bed_subset %>% filter(gene == target_gene,  tf == key_tfs) %>% as.data.frame
        all_var_in_tf <- all_var_in_tf[grepl(f_p('.*%s', key_features), all_var_in_tf$name),]
        transcript_data[key_features,1:10]
        if (length(setdiff(individual_snps_ids$snp_id, all_var_in_tf$snp_id)) != 0){
            ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
            
        }
        if (all(key_features_df[key_features_df$feature == key_features, c('feature_start', 'feature_end')] ==
                      transcript_data[key_features,c('feature_start', 'feature_end')]) == FALSE){
            ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
            
        }
        
        
        f_ASSERT( all(key_features_df[key_features_df$feature == key_features, c('feature_start', 'feature_end')] ==
                      transcript_data[key_features,c('feature_start', 'feature_end')]), "Features don't match")
        
        f_ASSERT(length(setdiff(individual_snps_ids$snp_id, all_var_in_tf$snp_id)) == 0, '')
        individual_snps_in_tf = individual_snps_ids$snp_id

        if (grepl('snpOnly', batch_name)){
            individual_snps_in_tf = intersect( individual_snps_in_tf, subset(all_var_in_tf, variant_type == 'common_snp' )$snp_id )
        }

                
        rownames(all_var_in_tf) = all_var_in_tf$snp_id
        recalculate_value = sum(all_var_in_tf[individual_snps_in_tf, 'impact'] * vcf_df[individual_snps_in_tf, 'add_score'])

        ##Compare the recalculate value and preprocessed value.
        if (recalculate_value != 0){
            digit_num = f_get_number_of_digits(transcript_data[key_features, individual_id])
            f_ASSERT( abs( round(recalculate_value, digits = digit_num) - transcript_data[key_features, individual_id] ) <= 1^-digit_num , f_p('Preprocess is not eaqual to recalculatio %s', individual_id))
            flog.info('Test past %s for %s variants in %s, %s are rare', recalculate_value, length(individual_snps_in_tf), individual_id, sum( all_var_in_tf[individual_snps_in_tf, 'variant_type'] == 'rare_var')  )
        }else{
            flog.info('recalculate_value is 0')
        }
    }

}


f_get_number_of_digits <- function(I){
    nchar(str_replace(as.character(I), '.*[.]', ''))
}


f_check_missing_genes <- function( merge_input, mode1, mode2){
    mode1_perf = merge_input %>% filter(mode == mode1)
    mode2_perf = merge_input %>% filter(mode == mode2)
    cat(f_p('%s - %s', mode1, mode2 ), '\n')
    setdiff(mode1_perf$gene, mode2_perf$gene)
    cat(f_p('%s - %s', mode2, mode1 ), '\n')
    setdiff(mode2_perf$gene, mode1_perf$gene)    
}

f_break_long_line_into_multiple <- function(long_str, wrap_len = 30){
    gsub(f_p('(.{1,%s})(\\s|$)', wrap_len), '\\1\n', long_str)
}


f_compare_modes_in_diff_batches <- function(batch_A, batch_B, mode_index1, mode_index2, mode_list, chr_list, perf_thres = 0.05, rsync_flag = F, performance_col = 'performance', add_batch_name = TRUE, debug = FALSE){

    if (f_judge_debug(debug)){
        ##:ess-bp-start::browser@nil:##
        browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
    }
    
    performance_merge1 = data.frame()
    performance_merge2 = data.frame()
    for (chr_str in chr_list){
        dir.create(f_p('./data/%s/rnaseq/%s/', batch_A, chr_str))
        return1 =f_summary_regression_results(batch_A, chr_str = chr_str, mode_name = mode_list[mode_index1], return_features = T, rsync_flag = rsync_flag)
        performance_merge1 = rbind(performance_merge1, return1$performance)

        dir.create(f_p('./data/%s/rnaseq/%s/', batch_B, chr_str))
        return2 =f_summary_regression_results(batch_B, chr_str = chr_str, mode_name = mode_list[mode_index2], return_features = T, rsync_flag = rsync_flag)
        performance_merge2 = rbind(performance_merge2, return2$performance)
    }

    if (add_batch_name){
        loc_batch_A = paste0(batch_A, ':',  mode_index1)
        loc_batch_B = paste0(batch_B, ':',  mode_index2)
    }else{
        loc_batch_A = mode_index1
        loc_batch_B = mode_index2
    }

    performance_merge1$mode = loc_batch_A
    performance_merge2$mode = loc_batch_B
    
    shared_cols=intersect(colnames(performance_merge1), colnames(performance_merge2))

    performance_merge = rbind(performance_merge1[, shared_cols], performance_merge2[, shared_cols])
    performance_merge$performance = performance_merge[,performance_col]
    
    dim(performance_merge2)
    dim(performance_merge1)

    setdiff(rownames(performance_merge1), rownames(performance_merge2) )
    setdiff(rownames(performance_merge2), rownames(performance_merge1) )

    table(performance_merge$mode)
    dot_plot=f_plot_performance_and_stats_test(subset(performance_merge, performance > perf_thres), loc_batch_A, loc_batch_B)
    compare_results = f_compare_improvment_for_two_groups(loc_batch_B, loc_batch_A, subset(performance_merge, performance > perf_thres), thres = 0.01, return_flag = F)

    #final_plot <-arrangeGrob( arrangeGrob(dot_plot, compare_results$diff, nrow =1),  compare_results$overall_density, compare_results$density, nrow = 3)
    final_plot <-arrangeGrob( arrangeGrob(dot_plot, compare_results$diff, nrow =1), compare_results$density, nrow = 2)
    plot(final_plot)

    return (performance_merge)
}



f_plot_values_against_performance <- function(plot_performance, target_col, test_stats = F){

    ##Check the performance against the mean expression values and variance of gene variance.    
        
    plot_performance$target = plot_performance[[target_col]]
    mean_cor = cor.test(plot_performance$performance, plot_performance$target)
    print(mean_cor)

    mean_cor$p.value = f_get_exact_pvalue(mean_cor)
    
    p <- ggplot(plot_performance, aes(performance, target)) + geom_point(alpha = 0.2, color ='#00BFC4') +
        annotate("text", label = f_p('P-value: %.2e\nCorr: %.2f', as.numeric(mean_cor$p.value), mean_cor$estimate),
                 x = 0.1, y = 1.6, hjust = 0) + ylab(target_col)+ geom_smooth()

    if(test_stats){
        p = list(pvalue = mean_cor$p.value, coef = mean_cor$estimate)
    }

    return (p)
}



f_adjust_lable_positions <- function(label_data, x_lab, y_lab, ...){
    library(FField)
    x.fact <- 100 / max(label_data[x_lab])
    y.fact <- 100 / max(label_data[y_lab])

    ## Repel points
    coords <- FFieldPtRep(coords = cbind(label_data[x_lab] * x.fact, 
                               label_data[y_lab] * y.fact),
                ...)
    head(coords)

    label_data$x.t <- coords$x / x.fact
    label_data$y.t <- coords$y / y.fact
    return (label_data)
}


f_overlap_tf_features_with_f5_enhancers <- function(feature_bed, enhancer_gene_filter, good_gene_list){
    feature_enhancer_bed = feature_bed[grepl('enhancer', feature_bed$name), ]
    dim(feature_enhancer_bed)
    head(enhancer_gene_filter)

    overlap_df = f_bedtools(feature_enhancer_bed, enhancer_gene_filter, paras = '-wao')
    head(overlap_df)
    overlap_df_sub = as.data.frame(overlap_df)[c('V1', 'V2', 'V3', 'V4', 'V8')]

    head(overlap_df_sub)
    colnames(overlap_df_sub) =c('chr','start', 'end', 'name', 'matched_gene' )
    overlap_df_gene  = overlap_df_sub %>% mutate( gene = str_replace(name, '[.].*', '') ) %>%
        filter( matched_gene == gene | matched_gene == '.', gene %in% good_gene_list) %>%
        mutate( f5_hit = matched_gene == gene )

    print(table(overlap_df_gene$gene == overlap_df_gene$matched_gene))
    return (overlap_df_gene)
}



f_test_one_tf_complete_in_gene_processed_data <- function(target_gene, batch_name, full_hic_loc, full_hic_contact, return_list_tf, debug = F){

    if (debug){
        ##:ess-bp-start::browser@nil:##
        browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
        
    }
    
    ##Read the gene data.
    test_gene <- GENE(data = data.frame(), gene_name = target_gene, chr_str = chr_str, batch_name = batch_name)
    test_gene$read_data()
    test_gene$subset_features_to_snp_contain_region('TF', debug = F)
    test_data_flag=test_gene$read_test_data(test_batch = '800samples', debug = FALSE)
    dim(test_gene$data)
    gene_pos_list=test_gene$hic_within_1Mb(return_pos= TRUE)
    transcript_data = test_gene$data
    rownames(transcript_data) = make.names(paste0(transcript_data$type, '-',transcript_data$feature), unique = T)


    ##Identify the first key features.
    key_features_df <- return_list_tf$features %>% separate(name, into =c('gene', 'feature'), sep ='[|]') %>% filter(gene == target_gene)
    key_features = grep('Intercept',key_features_df$feature,value = T, invert = T)
    if (length(key_features) != 1){
        key_features = key_features[1]
    }
    key_tfs = transcript_data[key_features,'feature']
    key_tfs


    ##Get all the promoter HiC ids of the gene.
    promoter_data = transcript_data %>% filter(type == 'promoter')
    enhancer_data = transcript_data %>% filter(type == 'enhancer') 
    promoter_hic_ids = unique(c(promoter_data$hic_fragment_id, enhancer_data$pair))
    promoter_hic_ids


    ##For the key features, use TF peak to overlap with all the promoter and associated regions
    full_interaction_pair <- full_hic_contact %>% filter(pair %in% promoter_hic_ids, score > 0.4) %>% select(pair, hic_fragment_id)
    full_interaction_hic <- full_hic_contact %>% filter(hic_fragment_id %in% promoter_hic_ids, score > 0.4) %>% select(pair, hic_fragment_id)
    dim(full_interaction_pair)
    dim(full_interaction_hic)
    head(full_interaction_pair)
    head(full_hic_loc)
    full_unique_hic_ids = unique(c(full_interaction_pair$hic_fragment_id, full_interaction_hic$pair))

    full_hic_ids_1M = full_unique_hic_ids[ full_unique_hic_ids >= gene_pos_list[1] & full_unique_hic_ids <= gene_pos_list[2] ]

    interact_hic_bed  <- full_hic_loc %>% filter(name %in% full_hic_ids_1M, chr == 'chr22' ) 


    #Get the TF bed regions which overlapped with at least one SNP
    ##tf_peak_file = list.files('./data/raw_data/tf/encode_peaks/processed/', pattern = f_p('.*%s.*narrowPeak', tolower(key_tfs)) , full.names = T)
    tf_var_file = list.files('./data/445samples_region/output/tf_variation/all/chr22/', pattern = f_p('%s.*txt', tolower(key_tfs)) , full.names = T)
    tf_bed_raw = read.table(tf_var_file, header = T)
    sample_cols = grep('(NA|HG)', colnames(tf_bed_raw), value = T)
    informative_rows = rowSums(tf_bed_raw[, sample_cols] != '.') >= 0.05 * 446
    tf_bed = tf_bed_raw[informative_rows,c(3,1,2,4)]



    #Get the TFs overlap with the gene's hic fragments
    overlapped_bed = f_bedtools(tf_bed, interact_hic_bed, fun = 'intersect', paras = '-wb')
    colnames(overlapped_bed) = c('chr', 'feature_start', 'feature_end', 'name', 'chr2', 'hic_start', 'hic_end', 'hic_fragment_id')

    unique_tf_hic=unique(as.data.frame(overlapped_bed)[,c('feature_start', 'feature_end', 'hic_fragment_id')])
    unique_tf_hic$name = with(unique_tf_hic, paste(feature_start, feature_end, hic_fragment_id, sep='_'))


    ###
    features_in_processed <- transcript_data %>% filter(feature == key_tfs) %>% select(feature, feature_start, feature_end ,type, cor ,hic_fragment_id, pair)
    unique_feature_regions = unique(features_in_processed[,c('feature_start', 'feature_end' ,'hic_fragment_id')])
    unique_feature_regions$name = with(unique_feature_regions, paste(feature_start, feature_end, hic_fragment_id, sep='_'))

    print('+++++++++++++++++++')
    ##This is possible, because not all the peaks will overlap with an variant.
    if (length(setequal(unique_feature_regions$name, unique_tf_hic$name))){
        
        flog.info('Test past %s, %s tf peaks', key_tfs, nrow(unique_feature_regions) )
    }else{
        flog.error('Test failed %s', key_tfs)
        print(setdiff(unique_feature_regions$name, unique_tf_hic$name))
        print(setdiff(unique_tf_hic$name, unique_feature_regions$name))
    }

}
