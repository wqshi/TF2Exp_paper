
from p_project_metadata import *

my.f_print_list(['===='])
print project_dir

####Subset the chr_vcf file into regulatory regions.
batch_name = '462samples'
batch_output_dir = f_get_batch_output_dir(batch_name)
#chr_num = '1'
kg_dir = '%s/data/raw_data/wgs/1kg/' % project_dir


#my.f_ensure_make_dir(output_dir)


#zcat ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | sed '/^[0-9]/s/^/chr/g'| head -n 10000 | bedtools intersect -header -a stdin -b encode-gm12878-mergeall.slop100Peak | bgzip > chr1.vcf.gz



###Subset each individual######12
sample_df = pd.read_csv('%s/output/sample.list' % batch_output_dir, sep = '\t', header = None)
print sample_df.head()
print sample_df.shape
sample_list = sample_df[1].tolist()

print 'NA12878' in sample_list

chr_list = []
for i in range(1,23):
    chr_list.append('%s' % i)
#for sample_id in sample_list:
#    print '===========%s============' % sample_id
#    #my.f_ensure_make_dir('%s/chr_vcf_files/%s/' % (kg_dir, sample_id))
#    for chr_num in chr_list:
#       bcf_subset_cmd = "~/packages/bcftools/bcftools-1.2/bcftools view -c1 -Ov -s %s %s | grep -v -i 'MULTI_ALLELIC\|ALU' | gzip > %s/chr_vcf_files/chr%s/%s.vcf.gz " % (sample_id, node_chr_vcf, output_dir, chr_num, sample_id)
#        my.f_shell_cmd(bcf_subset_cmd)
chr_list.append('X')
chr_list.append('Y')

chr_list = ['22']
for chr_num in chr_list:
    my.f_ensure_make_dir('%s/chr_vcf_files/chr%s/' % (batch_output_dir, chr_num))
    chr_vcf_file= "%s/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" % ( kg_dir, chr_num)
    subset_cmd = "zcat %s |  sed '/^[0-9]/s/^/chr/g' |bedtools intersect -header -a stdin -b %s/encode-gm12878-mergeall.slop100Peak | bgzip > %s/chr%s.vcf.gz" % (chr_vcf_file, kg_dir, kg_dir, chr_num)
    my.f_shell_cmd(subset_cmd)
    node_chr_vcf = '%s/chr%s.vcf.gz' % (kg_dir, chr_num)
    index_cmd = 'tabix -f -p vcf %s ' % node_chr_vcf
    my.f_shell_cmd(index_cmd)

quit()
def process_one_sample(sample_id, chr_list):
    #import ipdb; ipdb.set_trace()
    print '===========%s============' % sample_id
    for chr_num in chr_list:
        node_chr_vcf = '%s/chr%s.vcf.gz' % (kg_dir, chr_num)
        bcf_subset_cmd = "~/packages/bcftools/bcftools-1.2/bcftools view -c1 -Ov -s %s %s | grep -v -i 'MULTI_ALLELIC\|ALU' | awk -F'\t' '$1 ~ /^chr/{print $1,$2,$NF,$4,$5}' OFS='\t' | grep -v -P '\t<' | gzip > %s/chr_vcf_files/chr%s/%s.vcf.gz " % (sample_id, node_chr_vcf, batch_output_dir,chr_num, sample_id)
        my.f_shell_cmd(bcf_subset_cmd)

print sample_list

#Add NA12878 to the dataset.
#process_one_sample('NA12878', chr_list)
process_one_sample('HG00101', ['1'])

from joblib import Parallel, delayed  
import multiprocessing

# what are your inputs, and what operation do you want to 
# perform on each input. For example...
num_cores = 3 #multiprocessing.cpu_count()-4
print num_cores



for loc_chr in chr_list:
    break
    results = Parallel(n_jobs=num_cores)(delayed(process_one_sample)(sample_id, [loc_chr]) for sample_id in sample_list)






