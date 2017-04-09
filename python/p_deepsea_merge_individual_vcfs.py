
from p_project_metadata import *
import os
my.f_print_list(['===='])
print project_dir
import pandas as pd



####Subset the chr_vcf file into regulatory regions.
#batch_name = '462samples'
batch_name = '800samples'
batch_output_dir = f_get_batch_output_dir(batch_name)
#chr_num = '1'
print batch_output_dir

chr_num = '22'

def f_merge_vcf_files_for_one_chr(chr_num, batch_output_dir):


    chr_str    = 'chr%s' % chr_num
    print chr_str

    vcf_dir = '%s/chr_vcf_files/chr%s/' % (batch_output_dir, chr_num)

    raw_vcf_files = os.listdir(vcf_dir)

    individual_vcf_files = my.grep_list('(HG|NA)[0-9]*.vcf.gz', raw_vcf_files)

    #my.f_print_list(individual_vcf_files)
    
    vcf_gz_file = individual_vcf_files[1]

    vcf_list = []
    for vcf_gz_file in individual_vcf_files:
        try:
            vcf_data = pd.read_csv( '%s/%s' %(vcf_dir, vcf_gz_file), compression = 'gzip', sep = '\t', header = None)
        except:
            print 'Error in', vcf_gz_file
        vcf_list.append(vcf_data)


    merge_data = pd.concat(vcf_list)
    merge_data.columns = ['chr', 'start', 'genotype', 'ref', 'alt']

    print merge_data.shape

    reduced_data = merge_data.drop_duplicates(cols = ['chr', 'start', 'ref', 'alt'] )

    reduced_data.sort(columns=['chr', 'start'], inplace = True)

    print reduced_data.shape
    print vcf_dir
    reduced_data.to_csv('./data2/vcf/%s.vcf' % (chr_str), header = False, index = False, sep = '\t')


#for chr_num in range(1, 23):
f_merge_vcf_files_for_one_chr(chr_num, batch_output_dir)

from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()-2
print num_cores

#Parallel(n_jobs=6)(delayed(f_merge_vcf_files_for_one_chr)(chr_num, batch_output_dir) for chr_num in chr_list)



