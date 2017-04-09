
source('~/R/s_function.R', chdir = T)

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');
# base.dir = '.';


snp_file = './data/raw_data/wgs/1kg/additive_358samples/chr22.bed'
snp_bed = read.table(snp_file, header = T)
head(chr_data)



## Settings
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS




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
dim(snp_bed)
selected_samples = read.table('./data/raw_data/samples358.ped')
rm_cols=which(colnames(snp_bed) %in% selected_samples$V1 == FALSE)

snp_sample_order = colnames(snp_bed) [colnames(snp_bed) %in% selected_samples$V1]
head(snp_bed)
rownames(snp_bed) = make.names(snp_bed$SNP, unique = T)
snps = SlicedData$new();
snps$CreateFromMatrix(as.matrix(snp_bed[, snp_sample_order]))
#snps$fileDelimiter = "\t";      # the TAB character
#snps$fileOmitCharacters = "NA"; # denote missing values;
#snps$fileSkipRows = 1;          # one row of column labels
#snps$fileSkipColumns = rm_cols;       # one column of row labels
#snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
#snps$LoadFile(snp_file);


expression_file_name = './data/358samples_snyder_norm/rnaseq/GEUVADIS.Gene.DATA_MATRIX'
rna_seq_raw = read.csv(file = expression_file_name , sep=' ', header =TRUE)
head(rna_seq_raw)

## Load gene expression data



#gene$fileDelimiter = "\t";      # the TAB character
#gene$fileOmitCharacters = "NA"; # denote missing values;
#gene$fileSkipRows = 1;          # one row of column labels
#gene$fileSkipColumns = 1;       # one column of row labels
#gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
#gene$LoadFile(expression_file_name);

## Load covariates

## Run the analysis
snpspos = snp_bed[,c('SNP', 'chr', 'end')]
colnames(snpspos) = c('snpid', 'chr', 'pos')

library(stringr)
all.entrezgene = read.table('./data/raw_data/rnaseq/all.ensemble.genes.gene_start',sep='\t', header = TRUE, quote = "")
rownames(all.entrezgene) = all.entrezgene$ensembl_gene_id
all.entrezgene = subset(all.entrezgene, chromosome_name == '22')
head(all.entrezgene)

intersect_genes = intersect(str_replace( rna_seq_raw$gene, '[.].*', ''), all.entrezgene$ensembl_gene_id )
length(intersect_genes)
nrow(rna_seq_raw )

head(all.entrezgene)
my_genepos = all.entrezgene[ intersect_genes, c('ensembl_gene_id', 'chromosome_name', 'transcript_start', 'transcript_end')]
colnames(my_genepos) = c('geneid', 'chr', 'left', 'right')
my_genepos$chr = paste0('chr',  my_genepos$chr)
head(my_genepos)

show(gene)

show(snps)
show(gene)

head(snpspos)

gene = SlicedData$new();
rownames(rna_seq_raw) = str_replace(rna_seq_raw$gene, '[.].*', '')
rna_seq = rna_seq_raw[intersect_genes,snp_sample_order]
gene$CreateFromMatrix(as.matrix(rna_seq))
dim(rna_seq)

source('s_project_funcs.R')
show(snps)
show(gene)
head10(my_genepos)
head10(snpspos)

me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
                                        #cvrt = NULL,
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
cat('Detected local eQTLs:', '\n');
head(me$cis$eqtls)

input_eqtl = me$cis$eqtls

t_test_eqtls(input_eqtl, snpspos, my_genepos)
table(input_eqtl$gene)



t_test_eqtls <- function(input_eqtl, snpspos, my_genepos){
    random_rows = sample(rownames(input_eqtl), size = 100)

    library(RUnit)

    for (i in random_rows){
        cat(i, '\n')
        target_gene = input_eqtl[i, 'gene']
        target_snp = input_eqtl[i, 'snps']
        checkTrue( my_genepos[target_gene,'left'] - cisDist < snpspos[target_snp,'pos'], 'SNP > gene left' )
        checkTrue( my_genepos[target_gene,'right'] + cisDist > snpspos[target_snp,'pos'], 'SNP < gene right' )
        cat(target_gene,  my_genepos[target_gene,'left'] - cisDist, snpspos[target_snp,'pos'],  my_genepos[target_gene,'right'] + cisDist, '\n')
    }

}
## Plot the Q-Q plot of local and distant p-values

plot(me)

