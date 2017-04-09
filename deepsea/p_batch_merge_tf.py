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


#batch_name = '800samples'
batch_name = '462samples'
#chr_num_list =[22, 10, 15]
chr_num_list = ['X']

for chr_num in chr_num_list:
    cmd='python2.7 p_merge_tf_results.py --batch_name %s --chr_str chr%s --value_type diff' % (batch_name, chr_num)
    my.f_shell_cmd(cmd)
    cmd='python2.7 p_merge_tf_results.py --batch_name %s --chr_str chr%s --value_type ref' % (batch_name, chr_num)
    my.f_shell_cmd(cmd)


if my.f_get_server_name() == 'wqshi':
    if batch_name == '800samples':
        my.f_shell_cmd('scp $HOME/expression_var/data/%s/deep_result/all/chrMergeTF/*.gz shi@loire.cmmt.ubc.ca:/homed/home/shi/expression_var/data/800samples/deep_result/all/chr800/diff/' % (batch_name))
    else:
        my.f_shell_cmd('scp $HOME/expression_var/data/%s/deep_result/all/chrMergeTF/*.gz shi@loire.cmmt.ubc.ca:/homed/home/shi/expression_var/data/445samples_region/deep_result/all/chrMerge2/diff/' % (batch_name))












