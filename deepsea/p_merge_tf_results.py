from subprocess import *
from tempfile import mkdtemp
import sys
import os 
import os
home_dir = os.path.expanduser('~')
lib_dir = '%s/python/' % home_dir
import sys
sys.path.insert(0, lib_dir)
sys.path.insert(0, '%s/expression_var/python/' % home_dir)
import pandas as pd
import p_mymodule as my
from p_project_metadata import *

import argparse
parser = argparse.ArgumentParser(description='Extract the deepsea predictions of one tf.Deepsea prediction is sample based. The output of this script is TF based.')
if __doc__ is None:

    parser.add_argument('--batch_name', help = "batch_name", default =None)
    parser.add_argument('--chr_str', help = "chr", default =None)
    parser.add_argument('--value_type', help = "value type", default =None)
    args = parser.parse_args()
    
    batch_name = args.batch_name
    chr_str = args.chr_str
    value_type = args.value_type
else:
    
    batch_name = '445samples_region'
    chr_str = 'chr22'
    value_type = 'diff'

if my.f_get_server_name() == 'loire':
    batch_name = '445samples_region'
else:
    batch_name = args.batch_name

DEBUG=True

#vcf_file = '%s/deepsea/tests/data/%s.merge.head.vcf.gz'%(project_dir, chr_str)
vcf_file = '%s/data/%s/chr_vcf_files/chrMerge2/%s.vcf.gz'%(project_dir, batch_name, chr_str)
vcf_df = pd.io.parsers.read_csv(vcf_file, sep="\t", header=None, compression = 'gzip').ix[:,0:5]
print vcf_df.head()
vcf_df.columns = ['chr', 'pos', 'name', 'ref', 'alt']

print vcf_df.ix[:, 'ref'].head()
vcf_df.index = vcf_df.chr + '_' + vcf_df.pos.map(str) + '_' + vcf_df.ref + '_' + vcf_df.alt

deepsea_out = '%s/data/%s/deep_result/all/chrMergeTF/' % (project_dir, batch_name)

outdir = '%s/%s/' %(deepsea_out, chr_str)

logging.info('Out dir: %s', outdir)

print vcf_df.head()

peak_file_df_rmdup = f_get_peak_file_df_rmdup(project_dir, version = 'processed')

predictors = pd.read_csv('%s/deepsea/resources/predictor.names' % project_dir, header = None)
print predictors.head()
gm12878_predictors = f_add_suffix_on_duplicates(my.grep_list('gm12878', predictors.ix[:,0]))
logging.info('GM12878 features: %s' % len(gm12878_predictors))
for deepsea_col in gm12878_predictors:
    vcf_df[deepsea_col] = 0


target_tf = 'CTCF'
print filter(lambda x:re.search(r'GM12878[|]%s[|]'% target_tf, x, re.IGNORECASE), gm12878_predictors)

if f_judge_debug(DEBUG):
    import ipdb; ipdb.set_trace()

for loc_tf in peak_file_df_rmdup.tf:
    pred_file = '%s/%s.out.%s' %(outdir, loc_tf, value_type)
    target_tf = peak_file_df_rmdup.ix[loc_tf, 'deepsea_tf']
    print pred_file
    if os.path.isfile(pred_file):
        
        pred_data = pd.io.parsers.read_csv(pred_file, sep=",", header=0)
        pred_data.index = pred_data.chr + '_' + pred_data.pos.map(str) + '_' + pred_data.ref + '_' + pred_data.alt
        assert set(pred_data.index) < set(vcf_df.index), 'index error'
        assert len(set(pred_data.index)) == pred_data.shape[0], 'Duplicated index'
        #pred_data.columns = f_add_suffix_on_duplicates(pred_data.columns)

        for deepsea_col in my.grep_list('.*%s[|]' % target_tf, gm12878_predictors):
            vcf_df.ix[pred_data.index, deepsea_col] = pred_data.ix[:, deepsea_col]
    else:
        logging.info('Missing %s deepsea output' % loc_tf)

assert all(vcf_df.chr == chr_str), 'Error in chr'

print vcf_df.ix[:,1:10].head()
zero_variants=(vcf_df.ix[:, gm12878_predictors].sum(axis = 1) == 0).sum()
logging.info('%s out of %s are zero' % (zero_variants, vcf_df.shape[0]))


vcf_df.columns = [re.sub('None.[0-9]*', 'None', col_name) for col_name in vcf_df.columns]

print pd.isnull(vcf_df.pos).sum()

assert pd.isnull(vcf_df.pos).sum() == 0, 'No null positions'
print vcf_df.shape


print vcf_df.chr
vcf_df.index = range(vcf_df.shape[0])

outfile = '%s/data/%s/deep_result/all/chrMergeTF/%s.%s' % (project_dir, batch_name, chr_str, value_type)
vcf_df.to_csv( outfile, sep = ',',float_format='%.4e')
my.f_shell_cmd('gzip -f %s' % (outfile) )






