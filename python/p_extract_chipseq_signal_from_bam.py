
#Use homer to extract read density from bam files.
from p_project_metadata import *
import p_mymodule as my


if (my.f_get_server_name() == "loire"):
    head_dir="/homed/home/shi/anthony/tfbs_chipseq/ENCODE/dnase/"
    node_dir="/homed/home/shi/anthony/tfbs_chipseq/ENCODE/dnase/node_dir"
else:
    head_dir="/home/shi/projects/expression_var/data/raw_data/tf/embl_data"
    node_dir="/state/partition1/shi/tmp/"

tf_peak = {'PU1':'haib-gm12878-pu1.narrowPeak', 'RPB2':'haib-gm12878-pol2.narrowPeak', 'CTCF':'sydh-gm12878-ctcf.narrowPeak'}

#tf_list = ['RPB2', 'PU1']
tf_list = ['CTCF']


from joblib import Parallel, delayed
import multiprocessing

num_cores = 2 #multiprocessing.cpu_count()-4
print num_cores
#results = Parallel(n_jobs=num_cores)(delayed(process_one_sample)(sample_id, [new_chr_name]) for sample_id in sample_list)  


def f_process_one_CTCF(loc_bam, head_dir, node_base_dir):
    #import ipdb; ipdb.set_trace()
    individual_id = loc_bam.split('-')[1]+'_'+loc_bam.split('.')[0].split('-')[3]
    node_dir = node_base_dir + '/' + individual_id
    my.f_ensure_make_dir(node_dir)
    add_chr_cmd = "samtools view -H %s/%s | sed -e 's/SN:\([0-9XY]\)/SN:chr\\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - %s/%s > %s/%s" % (head_dir, loc_bam, head_dir, loc_bam, node_dir, loc_bam)
    my.f_shell_cmd(add_chr_cmd)
    individual_id = loc_bam.split('-')[1]+'_'+loc_bam.split('.')[0].split('-')[3]
    mkdir_cmd = 'makeTagDirectory %s/%s %s/%s' % (node_dir, individual_id, node_dir, loc_bam)
    my.f_shell_cmd(mkdir_cmd)
    copy_cmd = 'cp -r %s/%s %s; rm -r %s' % (node_dir, individual_id, head_dir, node_dir)
    my.f_shell_cmd(copy_cmd)

para_flag = True

for loc_tf in tf_list:

    data_dir = head_dir + '/' + loc_tf

    bam_list = my.f_grep_files_from_dir(data_dir, 'embl.*.bam', path=False)
    my.f_print_list(bam_list)
    loc_bam = bam_list[0]
    
    
    if para_flag == True:
        for loc_bam in bam_list:
            f_process_one_CTCF(loc_bam, data_dir, node_dir )
    else:
        Parallel(n_jobs=num_cores)(delayed(f_process_one_CTCF)(loc_bam, data_dir, node_dir) for loc_bam in bam_list)
    
    dir_list = my.f_grep_files_from_dir(data_dir, 'NA.*', path=True)
    my.f_print_list(dir_list)

    data_dirs=' '.join([loc_dir + '/' for loc_dir in dir_list])
    tf_dir = '%s/data/raw_data/tf/encode_peaks/processed/' % project_dir
    annotate_cmd = 'annotatePeaks.pl %s/%s hg19 -size given  -d %s -noann > %s/output.file' %(tf_dir, tf_peak[loc_tf], data_dirs, data_dir)
    my.f_shell_cmd(annotate_cmd)




#annotatePeaks.pl /homed/home/shi/expression_var/data/raw_data/tf/encode_peaks/processed/haib-gm12878-pu1.narrowPeak hg19 -size given  -d NA10851-PU1-Rep1/ NA10852-PU1-Rep1/ -noann > output.file

#RNA-seq


