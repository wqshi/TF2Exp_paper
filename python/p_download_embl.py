import os
import socket
import urllib
import re
import sys
import pandas as pd
import subprocess
import p_mymodule as my
import logging

server_name = socket.gethostname()
print server_name
logging.getLogger().setLevel(logging.DEBUG)
import p_mymodule as my




if (server_name == "loire"):
    head_dir="/homed/home/shi/anthony/tfbs_chipseq/ENCODE/dnase/"
   
    test_flag = True

else:
    head_dir="/home/shi/projects/chipseq_snp/data2/encode/"
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
        file_suffix=re.match(r"[a-zA-Z0-9_]*\.(.*)",file_name).group(1)
        #print file_suffix
    
        match_object=re.match(r".*(Rep[1-9]).*",file_name,flags=re.IGNORECASE)
        if match_object or len(file_names):
            if match_object:
                output_name=  prefix + "-"+match_object.group(1) + "." + file_suffix
            else:
                output_name=  prefix + "-Rep%s"%i + "." + file_suffix
                i = i + 1
        else:
            output_name = prefix + "." + file_suffix

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


embl_number = '1414'
index_url = 'http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-%s/E-MTAB-%s.sdrf.txt' % (embl_number, embl_number)
data_index_file = head_dir + 'embl_data.index'
urllib.urlretrieve(url=index_url , filename= data_index_file )
index_data = pd.read_csv(data_index_file, sep = '\t')

lab = 'embl'

logging.info('EMBL: %s' % embl_number)

if embl_number == '1884':

    tf_list = ['pu1', 'myc']
    cell_list = list(set(index_data['Source Name'].str.replace('NA', 'gm').tolist())) 
    
    my.f_print_list(cell_list)


    reload(my)


    for tf in tf_list:
        for cell in cell_list:
            data_dir = os.path.join(head_dir, cell)
            my.f_ensure_make_dir(data_dir)
            data_pattern = '%s.*%s' % (cell.replace('gm','na'), tf)
            fastq_prefix=my.f_create_file_prefix(cell, tf, lab, 'Rep1')
            fastq_field = 30
            #download_state=f_grep_wget_from_given_embl_file(data_index_file, data_pattern, data_dir, fastq_prefix, download_col = fastq_field, test_flag = test_flag, quiet = True, debug = False)
        
        
            bam_field = 36
            bam_prefix = my.f_create_file_prefix(cell, tf, lab, 'Rep2')
            #import ipdb; ipdb.set_trace()
            download_state=f_grep_wget_from_given_embl_file(data_index_file, data_pattern, data_dir, bam_prefix, download_col = bam_field, test_flag = test_flag, quiet = True, debug = False)
        
else:

    tf_list = ['cebpa','foxa1', 'hnf4a']
    #cell_list = list(set(index_data['Characteristics[strain]'].str.lower().tolist())) 
    cell_list = ['c57bl', 'aj', 'cast', 'spret']
    my.f_print_list(cell_list)

    
    for tf in tf_list:
        for cell in cell_list:
            data_dir = os.path.join(head_dir, cell)
            my.f_ensure_make_dir(data_dir)
            data_pattern = '%s.*%s' % (cell, tf)
            fastq_prefix=my.f_create_file_prefix(cell, tf, lab)
            download_state=f_grep_wget_from_given_embl_file(data_index_file, data_pattern, data_dir, fastq_prefix, download_pattern = 'ftp.*fastq.gz', test_flag = test_flag, quiet = True, debug = False)
        
        
#[Thu Apr  2 16:23:10 2015] p p_run_cluster_sep.py embl-download p_download_ebml.py 1 1 1 1 1 1 1 1 1
#[Thu Apr  2 16:23:58 2015] p p_run_cluster_sep.py embl-download p_download_ebml.py 1 1 1 1 1 1 1 1 1
#[Thu Apr  2 16:34:56 2015] p p_run_cluster_sep.py embl-download-shi p_download_ebml.py 1 1 1 1 1 1 1 1 1





#[Thu Apr  9 11:25:59 2015] p p_run_cluster_sep.py embl-download-shi p_download_embl.py 1 1 1 1 1 1 1 1 1
#[Wed Jun  3 10:54:51 2015] p p_run_cluster_sep.py embl-download-shi p_download_embl.py 1 1 1 1 1 1 1 1 1
