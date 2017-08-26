####################Using 462sample_log as base to fit the 462_peer model ####################################
import os
home_dir = os.path.expanduser('~')
lib_dir = '%s/python/' % home_dir
import sys
sys.path.insert(0, lib_dir)
sys.path.insert(0, '%s/expression_var/python/' % home_dir)
import pandas as pd
import p_mymodule as my
import argparse
parser = argparse.ArgumentParser(description='Run the cv.glmnet train on batch chrs')
import numpy as np

print '==========',__doc__

if __doc__ is None:
    parser.add_argument("--norm", help="The normalization method [original/norm/gtex]", default="norm")
    parser.add_argument('--batch_name', help = "462samples or 54samples", default = '462samples' )
    parser.add_argument('--chr_batch', help = "[chr22/3chrs/all]", default = 'chr22' )
    parser.add_argument('--modes', help = "1SNP, 2All, 3TF, 4 SNPinTF, 5 AlltfShuffle,  6TFShuffle, 7random, 8noInteract, 9fakeInteract, 10TFfilterMinor, 11TFsnpMatch, 12TFaddPenalty", default = '1,2,3' )
    parser.add_argument('--add_YRI', help = "Whether to include YRI population", default = 'TRUE' )
    parser.add_argument('--add_penalty', help = "Whether use hic interaction to guide feature selection", default = 'FALSE')
    parser.add_argument('--test', help ='Test one gene or all genes [TRUE/FALSE]', default = 'FALSE')
    parser.add_argument('--other_info', help = "[normCor/tradR2]", default = 'tradR2')
    parser.add_argument('--full_run', help = "run all or only good ones", default = 'True')
    
    args = parser.parse_args()
    add_YRI = args.add_YRI
    batch_name = args.batch_name
    norm_mode = args.norm
    chr_batch = args.chr_batch
    test_flag = args.test
    modes_index = np.array(args.modes.split(',')).astype(int).tolist()
    other_info = args.other_info
    add_penalty = args.add_penalty
    full_run = (args.full_run == 'True')
else:
    chr_batch = 'chr22'
    #batch_name='462samples_quantile_rmNA'
    #batch_name='445samples_sailfish'
    #batch_name='445samples_rareVar'
    batch_name = '445samples_snpOnly'
    norm_mode = 'norm'
    test_flag='TRUE' 
    modes_index = [1,2,3]
    other_info='keepZero'
    add_YRI = 'TRUE'
    add_penalty='FALSE'
    full_run = False   

if chr_batch == 'chr22':
    chr_list = [22]
elif chr_batch  == '3chrs':
    chr_list = [1, 22, 10, 15]
elif chr_batch  == '2chrs':
    chr_list = [10, 22]
elif chr_batch == 'rest':
    chr_list = list(set( range(1,23) + ['X'] ))  # - set([1,10,15,22]))
elif chr_batch == 'chrX':
    chr_list = ['X']
elif chr_batch == 'chr6':
    chr_list = [6]
else:
    chr_list = chr_list
# "1SNP, 2All, 3TF, 4 SNPinTF, 5 AlltfShuffle,  6TFShuffle, 7random, 8TFaddInteract, 9fakeInteract, 10TFfilterMinor, 11TFsnpMatch, 12TFaddPenalty"
#TFsnpMatch: 
mode_list=[ 'SNP', 'All', 'TF','SNPinTF', 'AlltfShuffle','TFShuffle',  'random', 'TFaddInteract', 'fakeInteract', 'TFfilterMinor', 'TFsnpMatch', 'TFaddPenalty']

gene='gene'
add_miRNA='FALSE'
add_TF_exp='FALSE'

add_TF_exp_only='FALSE'
add_predict_tf='FALSE'
add_TF_exp_only='FALSE'
#add_YRI='TRUE'
chr_str='chr22'
population='all'
TF_exp_type='fakeTF'
model='cv.glmnet'
add_gm12878='TRUE'

#other_info='rmdup'#??
#other_info='hic8'#Filter enhacer regions with hic score 0.8
#other_info='rmCor'
#other_info='lambda500'
#other_info = 'old2TFinter'

#####################################
#Input parameters

import re
print('Add YRI: %s', add_YRI)

sample_num = re.sub('samples.*', '', batch_name)
if norm_mode == 'norm':
    if 'Peer' in other_info:
        new_batch='%ssamples_peer' % sample_num #This one removes the population/gender, 27 hidden factors.
    elif 'GTex' in other_info:
        new_batch='%ssamples_gtex_norm' % sample_num #This one removes the population/gender, 27 hidden factors.
    else:
        new_batch='%ssamples_snyder_norm' % sample_num #This one removes the population/gender
else:
    new_batch='%ssamples_snyder_original' % sample_num


#chr_list=[10, 2, 22]
#chr_list=[22]


if norm_mode == 'norm':
    population = 'None'

for chr_num in chr_list:
    chr_str='chr%s' % chr_num
    print chr_str
    
    #mode_list=('All' 'SNP' 'SNPinTF' 'TF' 'AlltfShuffle' 'noInteract')
    #mode_list=('randomSNPinTF')
    #mode_list=('AlltfShuffle' 'AllsnpShuffle')
    #mode_list=['All', 'SNPinTF', 'random', 'AlltfShuffle']
    #mode_list=['AlltfShuffle']
    for new_batch_random in [mode_list[i-1] for i in modes_index]:
        
        run_cmd = 'sh s_start_cluster_gene_job.sh %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s; echo done' % (batch_name, test_flag, model, chr_str, gene, add_miRNA, add_TF_exp,  add_penalty, add_TF_exp_only, add_predict_tf, add_YRI, population, TF_exp_type, add_gm12878, new_batch, new_batch_random, other_info)
        my.f_shell_cmd(run_cmd)
        







