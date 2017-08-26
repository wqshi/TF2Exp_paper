
from p_project_metadata import *
#chr_list = [1..22] defined in the p_project_metadata
import argparse
reload(my)

parser = argparse.ArgumentParser(description='Extract the sample vcf from TF binding regions ')
print '==========',__doc__

if __doc__ is None:
    parser.add_argument("--chr_batch",     help="The chr wanted to compute[all/chr22]", default="chr22")
    parser.add_argument('--batch_name', help = "462samples or 54samples or 800samples", default = '462samples' )
    parser.add_argument('--test', help = "Run one[TRUE/FALSE]", default = 'FALSE' )
    args = parser.parse_args()
    if args.chr_batch != 'all':
        chr_list = ['X', 'Y']
    else:
        chr_list = chr_list
    batch_name = args.batch_name
    test_flag = args.test == 'TRUE'
else:
    chr_list = ['22']
    #batch_name = '462samples'
    batch_name = '800samples'
    test_flag = True

batch_output_dir = f_get_batch_output_dir(batch_name)
kg_dir = '%s/data/raw_data/wgs/1kg/' % project_dir

kg_samples = pd.read_csv('%s/integrated_call_samples_v3.20130502.ALL.panel' % kg_dir, sep='\t')

###Subset each individual######12
sample_df = pd.read_csv('%s/output/sample.list' % batch_output_dir, sep = '\t', header = None)
print sample_df.head()
print sample_df.shape
sample_list = sample_df[1].tolist()
sample_list=list(set(kg_samples.sample.tolist()).intersection(set(sample_list)))    
print 'Length %s' % len(sample_list)
for chr_num in chr_list:
    my.f_ensure_make_dir('%s/chr_vcf_files/chr%s/' % (batch_output_dir, chr_num))
    chr_vcf_file= "%s/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" % ( kg_dir, chr_num)

    if my.f_get_file_size('%s/chr%s.vcf.gz'%(kg_dir, chr_num)) > 1000: #> 1000k
        continue
    
    subset_cmd = "zcat %s |  sed '/^[0-9XY]/s/^/chr/g' |bedtools intersect -header -a stdin -b %s/encode-gm12878-mergeall.slop100Peak | bgzip > %s/chr%s.vcf.gz" % (chr_vcf_file, kg_dir, kg_dir, chr_num)
    my.f_shell_cmd(subset_cmd)
    node_chr_vcf = '%s/chr%s.vcf.gz' % (kg_dir, chr_num)
    index_cmd = 'tabix -f -p vcf %s ' % node_chr_vcf
    my.f_shell_cmd(index_cmd)

#sys.exit("Error message")

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
process_one_sample(sample_list[0], ['22'])

from joblib import Parallel, delayed  
import multiprocessing

# what are your inputs, and what operation do you want to 
# perform on each input. For example...
num_cores = 7 #multiprocessing.cpu_count()-4
print num_cores
print chr_list
print sample_list

if test_flag == False:
    for loc_chr in chr_list:
        results = Parallel(n_jobs=num_cores)(delayed(process_one_sample)(sample_id, [loc_chr]) for sample_id in sample_list)








