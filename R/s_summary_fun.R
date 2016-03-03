f_merge_data <- function(data_dir, file_pattern, quite = FALSE){

    feature_files = list.files(data_dir, pattern = file_pattern)

    cat('Total files', length(feature_files), '\n')
    feature_merge = data.frame()
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
        feature_merge = rbind(feature_merge, gene_feature)
    }

    return (feature_merge)
}


f_summary_regression_results <- function(batch_name, chr_str, mode_name, rsync_flag = TRUE){
    
    output_dir = f_p("./data/%s/rnaseq/%s/%s/", batch_name, chr_str, mode_name)

    #Rsync the features files back to the clustdell
    if (rsync_flag == TRUE){
        rsync_cmd = f_p("rsync -rav  --include '*enet*' --exclude '*' shi@clustdell.cmmt.ubc.ca:/home/shi/expression_var/data/%s/rnaseq/%s/%s/   ~/expression_var/data/%s/rnaseq/%s/%s/", batch_name, chr_str, mode_name, batch_name, chr_str, mode_name)
        system(rsync_cmd)
    }else{
        flog.info('Skip rsync...')
    }
    
    
    features_df=f_merge_data(output_dir, '.*enet.feature', quite = FALSE)
    
    performance_df = f_merge_data(output_dir, '.*enet$')
    head(performance_df)
    
    head(features_df)
    selected_feature_freq=as.data.frame(table(str_replace(features_df$name, '[|].*', '')))
    rownames(selected_feature_freq ) = selected_feature_freq$Var1
    performance_df$seleted_features = selected_feature_freq[performance_df$gene,'Freq']
    #hist(performance$selected_features)
    cat('Mean performance:',mean(performance_df$performance), '\n')
    #hist(performance_df$performance)

    good_predictions <- performance_df %>% filter(performance > 0.5)
    
    cat(nrow(good_predictions), 'out of', nrow(performance_df), 'have good predictions(Rsquared > 0.5)','\n')
    write.table(features_df, f_p('%s/features', output_dir), quote =FALSE, sep = '\t', row.names= FALSE)

    write.table(performance_df, f_p('%s/performance', output_dir), quote =FALSE, sep = '\t', row.names= FALSE)
    rownames(performance_df) = performance_df$gene
    return (performance_df)
}



f_collect_performance_in_mode_list <- function(batch_name, chr_str, mode_list, rsync_flag){
    collected_performance = NULL
    for (loc_mode in mode_list){
        cat('======', loc_mode, '========\n')
        performance = f_summary_regression_results( batch_name, chr_str, loc_mode, rsync_flag)
        if(is.null(collected_performance)){
            collected_performance = performance['performance']
            rownames(collected_performance) = performance$gene
            
        }else{
            shared_rows = intersect(performance$gene, rownames(collected_performance))
            length(shared_rows)
            collected_performance = cbind(collected_performance[shared_rows,], performance[shared_rows, 'performance'])
            rownames(collected_performance) = shared_rows
        }
    #head10(performance)
    
    }
    return (collected_performance)
}
