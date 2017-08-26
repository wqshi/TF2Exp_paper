####by Heather E. Wheeler 20150202####
##see runscripts/run_01_imputedDGN-WB_CV_elasticNet_chr*sh and qsub.txt for tarbell job submission scripts
source('s_project_funcs.R')
source('~/R/s_function.R', chdir = TRUE)
library(stringr)
source('s_gene_regression_fun.R')
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
args <- c('22', 0.5)
"%&%" = function(a,b) paste(a,b,sep="")

###############################################
### Directories & Variables


#snpset <- ".hapmapSnpsCEU" ##2015-02-02 results
#snpset <- ".wtcccGenotypedSNPs" ##2015-03-12 results
snpset <- "_1000G"

k <- 10 ### k-fold CV
tis <- "DGN-WB"  
chrom <- 22
chrname <- "chr" %&% chrom

##alpha = The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as
#(1-α)/2||β||_2^2+α||β||_1.
#alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.

alpha <- as.numeric(args[2]) #alpha to test in CV

################################################
### Functions & Libraries

library(glmnet)
#library(doMC) ##slower on tarbell than using 1 core, not sure why
#registerDoMC(10)
#getDoParWorkers()

################################################
source('s_gene_regression_fun.R')
all.entrezgene = f_get_all.entrezgene('./')

expdata=read.table(file = './data/462samples_quantile_rmNA/rnaseq/GEUVADIS.Gene.DATA_MATRIX', sep = ' ', header = T, quote = '')

row.names(expdata) = str_replace(expdata$gene, '[.].*', '')

head(all.entrezgene)
expdata$chr = all.entrezgene[row.names(expdata),'chromosome_name']
head10(expdata)
expdata = subset(expdata, chr == chrom )
dim(expdata)
expdata$chr=NULL
expdata$gene=NULL
t.expdata = t(expdata)


explist = colnames(t.expdata)

expsamplelist <- rownames(t.expdata) ###samples with exp data###
samplelist = expsamplelist
#bimfile <- gt.dir %&% "DGN.imputed_maf0.05_R20.8" %&% snpset %&% ".chr" %&% chrom %&% ".bim" ###get SNP position information###
#bim <- read.table(bimfile)
                
                        
exp.w.geno <- t.expdata[samplelist,] ###get expression of samples with genotypes###
explist <- colnames(exp.w.geno)

#gtfile <- gt.dir %&% "DGN.imputed_maf0.05_R20.8" %&% snpset %&% ".chr" %&% chrom %&% ".SNPxID"
#gtX <- scan(gtfile)
chr_str = 'chr22'
gtX = f_get_genetype_matrix(chr_str)

X <- gtX[,samplelist]
setdiff(samplelist, colnames(gtX))
bim = gtX[,c('CHR', 'SNP', 'POS')]
rownames(bim) = make.names(bim$SNP, unique = T)


X <- t(X) #transpose to match code below (one ID per row)


resultsarray = data.frame(gene='1',alpha=0,cvm=0,lambda.iteration=0,lambda.min=0,n.snps=0,R2=0,pval=0, caret_num=0, caret_R2=0, stringsAsFactors = F)
str(resultsarray)

en.dir = './'

weightcol = c("gene","SNP","refAllele","effectAllele","beta")
workingweight <- en.dir %&% tis %&% "_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_weights_chr" %&% chrom %&% "_" %&% date %&% ".txt"
write(weightcol,file=workingweight,ncol=5,sep="\t")

set.seed(1001) ##forgot to include in 2/2/15 run, should I re-run?
gencode = all.entrezgene
library(stringr)
i = 127




for(i in 1:length(explist)){
  cat(i,"/",length(explist),"\n")
  gene <- str_replace(explist[i], '[.].*', '')
  geneinfo <- all.entrezgene[gene,]
  c <- geneinfo$chromosome_name
  start <- geneinfo$gene_start - 1e6 ### 1Mb lower bound for cis-eQTLS
  end <- geneinfo$gene_end + 1e6 ### 1Mb upper bound for cis-eQTLs
  chrsnps <- subset(bim,bim[,1]==c) ### pull snps on same chr
  cissnps <- subset(chrsnps,chrsnps$POS >=start & chrsnps$POS <=end) ### pull cis-SNP info
 
  cisgenos <- ( X[,intersect(colnames(X),cissnps[,2])]) ### pull cis-SNP genotypes
 
  if(is.null(dim(cisgenos))){
    bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
  }else{
    minorsnps <- subset(colMeans(cisgenos), colMeans(cisgenos,na.rm=TRUE)>0) ###pull snps with at least 1 minor allele###
    minorsnps <- names(minorsnps)
    cisgenos <- cisgenos[,minorsnps]
    ##cisgenos <- scale(cisgenos, center=T, scale=T)
    ##cisgenos[is.na(cisgenos)] <- 0
    if(is.null(dim(cisgenos)) | dim(cisgenos)[2] == 0){###effectively skips genes with <2 cis-SNPs
      bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
    }else{

      exppheno <- exp.w.geno[,gene] ### pull expression data for gene
      exppheno <- scale(exppheno, center=T, scale=T)  ###need to scale for fastLmPure to work properly
      
      exppheno[is.na(exppheno)] <- 0
      rownames(exppheno) <- rownames(exp.w.geno)
      
      
      ##run Cross-Validation over alphalist
      sum(is.na(cisgenos))  
      cisgenos_df = as.data.frame(cisgenos)
      cisgenos_df$gene.RNASEQ = as.vector(exppheno)     
      try_results<-try(
          {fit <- cv.glmnet(cisgenos,exppheno,nfolds=k,alpha=alpha,keep=T,parallel=F) ##parallel=T is slower on tarbell, not sure why
            caret_fit  <- f_caret_regression_return_fit(my_train = cisgenos_df, target_col = 'gene.RNASEQ', learner_name = 'glmnet', tuneGrid = NULL)}
      )
      if (class(try_results) == "try-error"){
          flog.info('Error in %s', genename)
          next
      }
      
      fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm)) ##pull info to find best lambda
      best.lam <- fit.df[which.min(fit.df[,1]),] # needs to be min or max depending on cv measure (MSE min, AUC max, ...)
      cvm.best = best.lam[,1]
      lambda.best = best.lam[,2]
      nrow.best = best.lam[,3] ##position of best lambda in cv.glmnet output
      
      ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best]) # get betas from best lambda
      ret[ret == 0.0] <- NA
      bestbetas = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
      names(bestbetas) = rownames(ret)[which(!is.na(ret))]
      pred.mat <- fit$fit.preval[,nrow.best] # pull out predictions at best lambda

    }
  }
  genename <- gene
  if(length(bestbetas) > 0){
    res <- summary(lm(exppheno~pred.mat))
    
    rsq <- res$r.squared
    pval <- res$coef[2,4]

    
    resultsarray <- rbind(resultsarray,
                          c(genename, alpha, cvm.best, nrow.best, lambda.best, length(bestbetas), rsq, pval, caret_fit$selected_features, caret_fit$max_performance))

    ### output best shrunken betas for PrediXcan
    bestbetalist <- names(bestbetas)
    bestbetainfo <- bim[bestbetalist,]
    betatable<- (cbind(bestbetainfo,bestbetas))
    betafile<-cbind(genename,betatable$SNP, betatable$bestbetas) ##output "gene","SNP","refAllele","effectAllele","beta"
    #write.table(t(betafile),file=workingweight,ncolumns=5,append=T,sep="\t") # t() necessary for correct output from write() function
  }else{
    library(futile.logger)
    flog.info('Empty beta %s', genename)
    str(resultsarray)
    resultsarray <- rbind(resultsarray, c(genename, NA,NA,NA,NA,0,NA,NA, caret_fit$selected_features, caret_fit$max_performance))

  }
  #write(resultsarray[gene,],file=workingbest,ncolumns=8,append=T,sep="\t")
  resultscol <- c("gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval", 'caret_num', 'caret_R2')
  colnames(resultsarray) = resultscol
  write.table(resultsarray,file='predixcan.out',quote=F,row.names=F,sep="\t")

}


resultsarray=read.table('predixcan.out', header = TRUE)

mean(resultsarray$R2, na.rm = T)
mean(resultsarray$caret_R2)
sum(resultsarray$caret_R2>0.1, na.rm = T)/nrow(resultsarray)
dim(resultsarray)
