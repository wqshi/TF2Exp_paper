import re
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
import time
from p_project_metadata import *



import argparse

parser = argparse.ArgumentParser(description='Extract the deepsea predictions of one tf.Deepsea prediction is sample based. The output of this script is TF based.')

print '==========',__doc__

if __doc__ is None:
    parser.add_argument("--chr_str", help="The chr wanted to compute", default="chr22")
    parser.add_argument("--target_mode", help="TF/All", default="TF")
    parser.add_argument('--batch_name', help = "462samples or 54samples", default = '445samples_sailfish' )
    parser.add_argument('--other_info', help = "minor mode", default = 'tradR2' )
    parser.add_argument('--last_gene', help = "Last gene name to check", default = 'ENSG00000269103' )
    args = parser.parse_args()
    chr_str = args.chr_str
    target_mode = args.target_mode
    batch_name = args.batch_name
    other_info = args.other_info
    last_gene = args.last_gene
else:
    target_mode = 'TF'
    other_info = 'tradR2'
    chr_str = 'chr22'
    batch_name = '445samples_sailfish'
    last_gene = 'ENSG00000269103'



loc_dir = 'rm.histone_model.cv.glmnet_add.penalty_population.None_new.batch.445samples.snyder.norm_batch.mode.TFMODE_other.info.normCor'.replace('TFMODE', target_mode).replace('normCor', other_info)

full_result_dir = '%s/data/%s/rnaseq/%s/%s' % (project_dir, batch_name, chr_str, loc_dir)
print project_dir
print full_result_dir

import time
start_time = time.time()


if my.f_get_server_name() != 'loire' and 'clustdell' not in my.f_get_server_name(): 
    time_interval = 120
else:
    time_interval = 0.25



while True:
    gene_output = my.f_grep_files_from_dir(full_result_dir, '%s*.enet' % last_gene)
    
    if len(gene_output) == 0:
        time.sleep(time_interval)
    else:
        logging.info( 'Check the output %s' % gene_output)
    
        interval = start_time - os.path.getmtime(gene_output[0])
    
        if interval < 0:
            time.sleep(10*time_interval)
            my.f_shell_cmd('Rscript3 %s/R/r_summary_features_in_one_mode.R --batch_name %s --target_mode %s --chr_str %s ' %(project_dir, batch_name, loc_dir, chr_str ))
            break
        else:
            time.sleep(time_interval)


    if time.time() - start_time > 10*time_interval:
        my.f_shell_cmd('Rscript3 %s/R/r_summary_features_in_one_mode.R --batch_name %s --target_mode %s --chr_str %s ' %(project_dir, batch_name, loc_dir, chr_str ))
        break


