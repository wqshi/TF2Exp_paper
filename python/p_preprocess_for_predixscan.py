
from p_project_metadata import *
kg_dir = '%s/data/raw_data/wgs/1kg/' % project_dir
import pandas as pd
import pybedtools


kg_dir = '%s/data/raw_data/wgs/1kg/' % project_dir
###Merge the Gene regions#####
#all_entrezgene = pd.read_csv('./data/raw_data/rnaseq/all.ensemble.genes.gene_start',sep='\t')
#all_entrezgene['chr'] = 'chr' + all_entrezgene['chromosome_name']
#bed_data = all_entrezgene.ix[:,['chr', 'transcript_start', 'transcript_end', 'ensembl_gene_id']]
#bed_obj = my.f_pd_to_bed(bed_data)
#merged_bed = bed_obj.sort().slop(g=pybedtools.chromsizes('hg19'), b= 1000000).merge(c=1, o='count')
#merged_bed.saveas('%s/gene.1MB.merge.bed' % kg_dir)



##Filter the vcf SNPs####

chr_list=['22', '21', 'X']

chr_vcf_file = '%s/ALL.head.vcf.gz' % kg_dir
chr_num = 'Test'

additive_dir = '%s/additive_445samples/' % kg_dir
my.f_ensure_make_dir(additive_dir)

for chr_num in chr_list:
    
    bcftools_dir = '~/packages/bcftools/bcftools-1.2/'
    chr_vcf_file= "%s/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" % ( kg_dir, chr_num)
    #my.f_shell_cmd(subset_cmd)
    #sample_file = './data/raw_data/samples370.ped'
    loc_sample_file = './data/raw_data/samples445.ped'
    sample_file = '%s/data/raw_data/samplesAll.ped' % project_dir
    sample_vcf_file = '%s/chr%s.vcf.gz' % (additive_dir, chr_num)

    if chr_num == 'Y':
        force_flag = '--force-samples'
    else:
        force_flag = ''
    
    #Use the MAF from local population
    loc_maf_cmd = "%s/bcftools view -c1 -Ov %s --samples-file %s %s | %s/bcftools filter -e'MAF<0.05' - | grep -v -i 'MULTI_ALLELIC\|ALU\|#' | awk '{print \"chr\"$1\"\t\"$2-1\"\t\"$2\"\t\"$3}' > %s/chr%s.vcf.loc" % (bcftools_dir, force_flag, loc_sample_file, chr_vcf_file, bcftools_dir, additive_dir, chr_num)
    print loc_maf_cmd
    my.f_shell_cmd(loc_maf_cmd)
    
    #Use the MAF from global population, other wise change of samples would affect the SNPs selected.

    
    bcf_subset_cmd = " %s/bcftools view -c1 -Ov %s --samples-file %s %s | awk '{if($0 !~ /^#/) print \"chr\"$0; else print $0}' | bedtools intersect -header -b %s/chr%s.vcf.loc -a stdin -wa | bgzip > %s " % ( bcftools_dir, force_flag, sample_file, chr_vcf_file, additive_dir, chr_num, sample_vcf_file)

    print bcf_subset_cmd
    my.f_shell_cmd(bcf_subset_cmd)

    plink_cmd = '%s/plink --vcf %s --recode A-transpose  --out %s/chr%s --noweb' % (kg_dir, sample_vcf_file, additive_dir, chr_num )
    my.f_shell_cmd(plink_cmd)

    format_cmd = "awk '{print \"chr\"$0\"\t\"$4-1}' %s/chr%s.traw | sed '1s/POS/end/' | sed '1s/-1/start/' | sed '1s/CHR/chr/'  | sed '1s/chrchr/chr/g'> %s/chr%s.bed" % (additive_dir, chr_num, additive_dir, chr_num)
    print(format_cmd)
    my.f_shell_cmd(format_cmd)
