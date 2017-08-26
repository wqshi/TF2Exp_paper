import os
#dir_path = os.path.dirname(os.path.realpath(__file__))
#print dir_path

import socket
import urllib
import re
import sys
import pandas as pd
import subprocess

home_dir = os.path.expanduser('~')
lib_dir = '%s/python/' % home_dir
import sys
sys.path.insert(0, lib_dir)

import p_mymodule as my
import logging

server_name = socket.gethostname()
print server_name
logging.getLogger().setLevel(logging.DEBUG)
import p_mymodule as my




if (server_name == "loire"):
    head_dir="/homed/home/shi/anthony/tfbs_chipseq/ENCODE/dnase/"
    test_flag = False
else:
    head_dir="/home/shi/projects/expression_var/data/raw_data/tf/embl_data"
    test_flag = False


def f_grep_wget_from_given_embl_file(index_file, pattern,  output_dir, prefix, download_pattern ,test_flag=False, quiet=False, debug = False):
    if debug == True:
        import ipdb; ipdb.set_trace()
        
    import urllib
    matched_lines=my.grep_file(pattern, index_file)
    if matched_lines == None:
        if quiet == False:
            print "-----------------------Warning--------------------------\nNo matching for the pattern %s in %s\n"%(pattern, index_file)
        return "failed"
    
    file_names=[ my.grep_list( download_pattern,  re.split('\t',line) )[0] for line in matched_lines]
    #print file_names

    i=1
    for file_name in file_names:

        data_url = os.path.dirname(file_name)
        file_name = os.path.basename(file_name)
        #file_suffix=re.match(r"[a-zA-Z0-9_]*\.(.*)",file_name).group(1)
        #print file_suffix

        tmp, file_suffix = os.path.splitext(file_name)
        
        match_object=re.match(r".*(Rep[1-9]).*",file_name,flags=re.IGNORECASE)
        if match_object or len(file_names):
            if match_object:
                output_name=  prefix + "-"+match_object.group(1) + file_suffix
            else:
                output_name=  prefix + "-Rep%s"%i  + file_suffix
                i = i + 1
        else:
            output_name = prefix + file_suffix

        output_file = output_dir + "/" + output_name

        if test_flag == False:
            urllib.urlretrieve(url=data_url + "/" + file_name, filename= output_file )
            
        print "Downlaod " + file_name + " " + data_url + ' ' + output_file
        
        match_object=re.match(r".*\Peak.gz$",file_name)
        if match_object:
            if test_flag == False:
                f_unzip_targz(output_file)
            print "Unzip " + output_name
        
    return "success"


embl_number = '3657'


#for embl_number in ['3656', '3657']:
index_url = 'http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-%s/E-MTAB-%s.sdrf.txt' % (embl_number, embl_number)
data_index_file = head_dir + 'embl_data.index'
urllib.urlretrieve(url=index_url , filename= data_index_file )
index_data = pd.read_csv(data_index_file, sep = '\t')

lab = 'embl'

logging.info('EMBL: %s' % embl_number)

print index_data.head()

if embl_number == '3657':
    tf_list = ['PU1', 'RPB2']
elif embl_number == '3656':
    tf_list = ['RNA']
else:
    print 'Wrong embl number'



if embl_number != '': #TF binding data

    
    cell_list = list(set(index_data['Characteristics[coriell id]'].str.replace('NA', 'NA').tolist())) 
    
    my.f_print_list(cell_list)
    for tf in tf_list:
        data_dir = os.path.join(head_dir, tf)
        my.f_ensure_make_dir(data_dir)
        for cell in cell_list:
            data_pattern = '%s_%s' % (cell.replace('NA',''), tf)
            #fastq_prefix=my.f_create_file_prefix(cell, tf, lab, 'Rep1')
            #fastq_field = 30
            #download_state=f_grep_wget_from_given_embl_file(data_index_file, data_pattern, data_dir, fastq_prefix, download_col = fastq_field, test_flag = test_flag, quiet = True, debug = False)
            bam_field = 36
            bam_prefix = my.f_create_file_prefix(cell, tf, lab)
            #import ipdb; ipdb.set_trace()
            download_state=f_grep_wget_from_given_embl_file(data_index_file, data_pattern, data_dir, bam_prefix, download_pattern = 'ftp.*bam', test_flag = test_flag, quiet = True, debug = False)
            if my.f_get_server_name() == 'loire':
                break
        
#[Fri Feb 10 14:05:19 2017] p p_run_cluster_sep.py download-shi p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 14:54:14 2017] p p_run_cluster_sep.py download-shi p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 14:55:41 2017] p p_run_cluster_sep.py download-shi p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 14:59:58 2017] p p_run_cluster_sep.py download-shi p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 15:01:41 2017] p p_run_cluster_sep.py download-shi p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 15:10:36 2017] p p_run_cluster_sep.py download-shi p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 15:12:28 2017] p p_run_cluster_sep.py download-shi p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 15:13:09 2017] p p_run_cluster_sep.py download-shi p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 15:55:08 2017] p p_run_cluster_sep.py download-shi p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 15:56:14 2017] p p_run_cluster_sep.py download-shi-tf p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 16:43:01 2017] p p_run_cluster_sep.py download-shi-tf p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 16:43:46 2017] p p_run_cluster_sep.py download-shi-rna p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
#[Fri Feb 10 20:49:00 2017] p p_run_cluster_sep.py download-shi-tf p_download_embl.py 1 1 1 1 1 1 1 1 1 1 1
