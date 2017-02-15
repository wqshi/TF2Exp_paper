
library(futile.logger)
library(stringr)
library(dplyr)
library(futile.logger)
head10 <- function(loc_data){

    row_num = min(nrow(loc_data), 10)
    col_num   = min(ncol(loc_data), 10)
    print(loc_data[1:row_num, 1:col_num])
}


f_heatmap_genes <- function(loc_data, output_file, seed_number = NULL){

    
    library(gplots)
    tiff(output_file, width = 7, height = 7, units = 'in', res = 300)
    if (!is.null(seed_number)){
        set.seed(seed_number)
    }
    subset_index = sample(rownames(loc_data), size = 2000)
    library('RColorBrewer')
    RdBu_color <- colorRampPalette(brewer.pal(n = 10, "RdBu"))
    heatmap.2(loc_data[subset_index,], trace = 'none', main ='Random 2000 genes', col = RdBu_color(256) )
    dev.off()
}

options(stringsAsFactors = FALSE)

f_check_na_rows <- function(input_data){
    na_index=which(rowSums(is.na(input_data)) > 0)
    flog.info('%s NA rows detected', length(na_index))
    return (na_index)
}


f_filter_NA_rm_zeroVar <- function(rna_seq_raw){
    rna_seq_raw_filter=rna_seq_raw[complete.cases(rna_seq_raw),]
    flog.info('Filter NA: %s out of %s', nrow(rna_seq_raw_filter), nrow(rna_seq_raw))

    f_check_na_rows(rna_seq_raw_filter)
    #rna_seq_raw_rmZero = rna_seq_raw_filter[rowSums(rna_seq_raw_filter == 0) < 0.25*ncol(rna_seq_raw),]
    near_zero = caret::nearZeroVar(t(rna_seq_raw_filter))
    if(length(near_zero) == 0){
        rna_seq_raw_rmZero = rna_seq_raw_filter
    }else{
        rna_seq_raw_rmZero = rna_seq_raw_filter[-(near_zero),]
    }

    flog.info('Filter Zero: %s out of %s', nrow(rna_seq_raw_rmZero), nrow(rna_seq_raw_filter))

    return (rna_seq_raw_rmZero)
}

f_find_NA_cols <- function(input_data){
    a=colSums(is.na(input_data))
    return (names(a)[a>1])
}



f_duplicated_cols <- function(df){
    return (duplicated(df) | duplicated(df[nrow(df):1,])[nrow(df):1])

}

f_judge_debug<-function(debug_flag){
    return (f_get_server_name() == 'loire' & debug_flag)
}



f_break_long_line_into_multiple <- function(long_str, wrap_len = 40){
    gsub(f_p('(.{1,%s})(\\s|$)', wrap_len), '\\1\n', long_str)
}

f_get_exact_pvalue <- function(A){
    2* pt(A$statistic,  df = A$parameter, lower.tail=FALSE)
}


f_row_scale <- function(df){
    return (t(scale(t(df))))
}

f_quantile <- function(input_data){

    data_quantile = normalize.quantiles(as.matrix(input_data))
    colnames(data_quantile) = colnames(input_data)
    rownames(data_quantile) = rownames(input_data)
    return (data_quantile)
}
