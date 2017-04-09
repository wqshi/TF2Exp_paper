
from subprocess import *
from tempfile import mkdtemp
import sys
import os 
home_dir = os.path.expanduser('~')
lib_dir = '%s/python/' % home_dir
import sys
sys.path.insert(0, lib_dir)
sys.path.insert(0, '%s/expression_var/python/' % home_dir)
import pandas as pd
import p_mymodule as my
from p_project_metadata import *


chr_list = [10, 15, 22]
output_dir = '%s/data/462samples/deep_result/all/chrMergeTF' % project_dir

for chr_num in chr_list:
    chr_str = 'chr' + str(chr_num)
    loc_output_dir = '%s/%s' %(output_dir, chr_str)
    loc_vcf_file  = '%s/data/462samples/chr_vcf_files/chrMerge2/%s.vcf.gz' % (project_dir, chr_str)
    cmd = 'python2.7 p_rundeepsea.py --vcf_file %s --out_dir %s' % (loc_vcf_file, loc_output_dir)
    my.f_shell_cmd(cmd)
    
    
