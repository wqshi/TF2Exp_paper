
source('~/R/s_function.R', chdir = T)
library(MatrixEQTL)
source('s_project_funcs.R')
source('~/R/s_ggplot2_theme.R', chdir = TRUE) 
f_call_eQTL <- function(snp_file, snp_test_file, selected_samples, expression_file_name, snpspos_file){
    useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
    snp_bed = read.table(snp_test_file, header = T)
                                        # Output file name
    output_file_name_cis = tempfile();
    output_file_name_tra = tempfile();

                                        # Only associations significant at this level will be saved
    pvOutputThreshold_cis = 1e-2;
    pvOutputThreshold_tra = 0;

                                        # Error covariance matrix
                                        # Set to numeric() for identity.
    errorCovariance = numeric();
                                        # errorCovariance = read.table("Sample_Data/errorCovariance.txt");

                                        # Distance for local gene-SNP pairs
    cisDist = 1e6;

    ## Load genotype data
    keep_cols=which(colnames(snp_bed) %in% selected_samples$V1)
    snp_sample_order = colnames(snp_bed) [colnames(snp_bed) %in% selected_samples$V1]

    snps = SlicedData$new();
                                        #snps$CreateFromMatrix(as.matrix(snp_bed[, snp_sample_order]))
    snps$fileDelimiter = "\t";      # the TAB character
    snps$fileOmitCharacters = "NA"; # denote missing values;
    snps$fileSkipRows = 1;          # one row of column labels
    snps$fileSkipColumns = 1;       # one column of row labels
    snps$fileSliceSize = 10000;      # read file in slices of 2,000 rows
    snps$LoadFile(snp_file);
    snps$ColumnSubsample(keep_cols - 1)

    

    ## Read SNP positions
    snpspos = read.table(snpspos_file, header = T)
    colnames(snpspos) = c('chr', 'snp', 'pos')
    snpspos = snpspos[,c('snp', 'chr', 'pos')]
    snpspos$pos = as.numeric(snpspos$pos)
    snpspos = snpspos[complete.cases(snpspos),]



    ##Read gene locations
    library(stringr)
    all.entrezgene = read.table('./data/raw_data/rnaseq/all.ensemble.genes.gene_start',sep='\t', header = TRUE, quote = "")
    rownames(all.entrezgene) = all.entrezgene$ensembl_gene_id
    head(all.entrezgene)


    rna_seq_raw = read.csv(file = expression_file_name , sep=' ', header =TRUE)
    intersect_genes = intersect(str_replace( rna_seq_raw$gene, '[.].*', ''), all.entrezgene$ensembl_gene_id )
    flog.info('Length %s', length(intersect_genes))

    my_genepos = all.entrezgene[ intersect_genes, c('ensembl_gene_id', 'chromosome_name', 'transcript_start', 'transcript_end')]
    colnames(my_genepos) = c('geneid', 'chr', 'left', 'right')
    my_genepos$chr = paste0('chr',  my_genepos$chr)



                                        #The gene object
    gene = SlicedData$new();
    rownames(rna_seq_raw) = str_replace(rna_seq_raw$gene, '[.].*', '')
    rna_seq = rna_seq_raw[intersect_genes,snp_sample_order]
    gene$CreateFromMatrix(as.matrix(rna_seq))
    dim(rna_seq)

    me = Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        ##cvrt = NULL,
        output_file_name     = output_file_name_tra,
        pvOutputThreshold     = pvOutputThreshold_tra,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos,
        genepos = my_genepos,
        cisDist = cisDist,
        pvalue.hist = "qqplot",
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);

    unlink(output_file_name_tra);
    unlink(output_file_name_cis);

    ## Results:

    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    cat('Detected local eQTLs:', length(me$cis$eqtls), '\n');

    return (list(eqtl=me$cis$eqtls, genepos = my_genepos, snpspos = snpspos))

}

t_test_eqtls <- function(){

    snp_file = './data/raw_data/wgs/1kg/additive_358samples/merge.snps'
    snp_test_file = './data/raw_data/wgs/1kg/additive_358samples/merge.snps.test'
    cisDist = 1e6
    selected_samples = read.table('./data/raw_data/samples358.ped')
    expression_file_name = './data/358samples_snyder_norm/rnaseq/GEUVADIS.Gene.DATA_MATRIX'
    snpspos_file = './data/raw_data/wgs/1kg/additive_358samples/snp_pos'

    
    return_list = f_call_eQTL(snp_test_file, snp_test_file, selected_samples, expression_file_name, snpspos_file)
    input_eqtl = return_list$eqtl
    my_genepos = return_list$genepos
    snpspos = return_list$snpspos
    
    random_rows = sample(rownames(input_eqtl), size = min(100, nrow(input_eqtl)))

    library(RUnit)
    rownames(snpspos) = make.names(snpspos$snp, unique = T)
    for (i in random_rows){
        target_gene = input_eqtl[i, 'gene']
        target_snp = input_eqtl[i, 'snps']
        checkTrue( my_genepos[target_gene,'left'] - cisDist < snpspos[target_snp,'pos'], 'SNP > gene left' )
        checkTrue( my_genepos[target_gene,'right'] + cisDist > snpspos[target_snp,'pos'], 'SNP < gene right' )
        cat(target_gene,  my_genepos[target_gene,'left'] - cisDist, snpspos[target_snp,'pos'],  my_genepos[target_gene,'right'] + cisDist, '\n')
    }

}


#t_test_eqtls()



f_find_best_peer_number <- function(snp_file, snp_test_file, selected_samples, expression_file_name, snpspos_file, peer_numbers, debug = FALSE){

    if (debug == TRUE){
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]))##:ess-bp-end:##
        
    }

    
    useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
    snp_bed = read.table(snp_test_file, header = T)
                                        # Output file name
    output_file_name_cis = tempfile();
    output_file_name_tra = tempfile();

                                        # Only associations significant at this level will be saved
    pvOutputThreshold_cis = 1e-2;
    pvOutputThreshold_tra = 0;

                                        # Error covariance matrix
                                        # Set to numeric() for identity.
    errorCovariance = numeric();
                                        # errorCovariance = read.table("Sample_Data/errorCovariance.txt");

                                        # Distance for local gene-SNP pairs
    cisDist = 1e6;

    ## Load genotype data
    keep_cols=which(colnames(snp_bed) %in% selected_samples$V1)
    snp_sample_order = colnames(snp_bed) [colnames(snp_bed) %in% selected_samples$V1]

    snps = SlicedData$new();
                                        #snps$CreateFromMatrix(as.matrix(snp_bed[, snp_sample_order]))
    snps$fileDelimiter = "\t";      # the TAB character
    snps$fileOmitCharacters = "NA"; # denote missing values;
    snps$fileSkipRows = 1;          # one row of column labels
    snps$fileSkipColumns = 1;       # one column of row labels
    snps$fileSliceSize = 10000;      # read file in slices of 2,000 rows
    snps$LoadFile(snp_file);
    snps$ColumnSubsample(keep_cols - 1)

    

    ## Read SNP positions
    snpspos = read.table(snpspos_file, header = T)
    colnames(snpspos) = c('chr', 'snp', 'pos')
    snpspos = snpspos[,c('snp', 'chr', 'pos')]
    snpspos$pos = as.numeric(snpspos$pos)
    snpspos = snpspos[complete.cases(snpspos),]



    ##Read gene locations
    library(stringr)
    all.entrezgene = read.table('./data/raw_data/rnaseq/all.ensemble.genes.gene_start',sep='\t', header = TRUE, quote = "")
    rownames(all.entrezgene) = all.entrezgene$ensembl_gene_id
    head(all.entrezgene)


    
    counts = rep(0, length(peer_numbers))
    i = 1
    for (peer_number in peer_numbers){

        flog.info('PEER factor number %s', peer_number)
        rna_seq_raw = read.csv(file = paste0(expression_file_name, '.peer', peer_number), sep=' ', header =TRUE)
        intersect_genes = intersect(str_replace( rna_seq_raw$gene, '[.].*', ''), all.entrezgene$ensembl_gene_id )
        flog.info('Length %s', length(intersect_genes))

        my_genepos = all.entrezgene[ intersect_genes, c('ensembl_gene_id', 'chromosome_name', 'transcript_start', 'transcript_end')]
        colnames(my_genepos) = c('geneid', 'chr', 'left', 'right')
        my_genepos$chr = paste0('chr',  my_genepos$chr)

        

        ##The gene object
        gene = SlicedData$new();
        rownames(rna_seq_raw) = str_replace(rna_seq_raw$gene, '[.].*', '')
        rna_seq = rna_seq_raw[intersect_genes,snp_sample_order]
        gene$CreateFromMatrix(as.matrix(rna_seq))
        dim(rna_seq)

        me = Matrix_eQTL_main(
            snps = snps,
            gene = gene,
            ##cvrt = NULL,
            output_file_name     = output_file_name_tra,
            pvOutputThreshold     = pvOutputThreshold_tra,
            useModel = useModel,
            errorCovariance = errorCovariance,
            verbose = TRUE,
            output_file_name.cis = output_file_name_cis,
            pvOutputThreshold.cis = pvOutputThreshold_cis,
            snpspos = snpspos,
            genepos = my_genepos,
            cisDist = cisDist,
            pvalue.hist = "qqplot",
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = FALSE);

        unlink(output_file_name_tra);
        unlink(output_file_name_cis);

                
        counts[i] = length(  unique( subset(me$cis$eqtls, pvalue < 1e-8)$gene ))
        i = i + 1
        cat('Analysis done found ', counts[peer_number], ' genes and ',  sum(me$cis$eqtls$FDR < 0.05), ' eQTLs in: ', me$time.in.sec, ' seconds', '\n')
    }
    return (counts)

}

t_test_peer_number <- function(){
    snp_file = './data/raw_data/wgs/1kg/additive_358samples/merge.snps'
    snp_test_file = './data/raw_data/wgs/1kg/additive_358samples/merge.snps.test'
    cisDist = 1e6
    selected_samples = read.table('./data/raw_data/samples358.ped')
    expression_file_name = './data/358samples_snyder_norm/rnaseq/GEUVADIS.Gene.DATA_MATRIX'
    snpspos_file = './data/raw_data/wgs/1kg/additive_358samples/snp_pos'

    
    
    peer_numbers = seq(from = 2, to = 5, 2)    
    return_list = f_find_best_peer_number(snp_test_file, snp_test_file, selected_samples, expression_file_name, snpspos_file, peer_numbers = peer_numbers, debug = F)

    peer_numbers = seq(from = 2, to = 50, 2)
    return_list = f_find_best_peer_number(snp_file, snp_test_file, selected_samples, expression_file_name, snpspos_file, peer_numbers = peer_numbers, debug = F)

    results_df = data.frame(PEER_number = peer_numbers, gene_num = return_list)

    ggplot(results_df, aes(PEER_number, gene_num)) + geom_point(color = 'blue') + xlab('Number of hidden factors in PEER') +
        ylab('CIS QTL genes') + theme_Publication(base_size = 14) + geom_line(color = 'blue')
    ggsave(f_p('./result/figures/peer_hidden_factors_eqtl.tiff'))
}


t_test_peer_number()





