
#Use homer to extract read density from bam files.
from p_project_metadata import *
import p_mymodule as my


if (my.f_get_server_name() == "loire"):
    head_dir="/homed/home/shi/anthony/tfbs_chipseq/ENCODE/dnase/"
else:
    head_dir="/home/shi/projects/expression_var/data/raw_data/tf/embl_data"



tf_list = ['RPB2', 'PU1']

loc_tf = 'RBP2'



tf_peak = {'PU1':'haib-gm12878-pu1.narrowPeak', 'RPB2':'haib-gm12878-pol2.narrowPeak'}

for loc_tf in tf_list:

    data_dir = head_dir + '/' + loc_tf

    bam_list = my.f_grep_files_from_dir(data_dir, '.*bam', path=False)
    my.f_print_list(bam_list)

    loc_bam = bam_list[0]
    for loc_bam in bam_list:
        individual_id = loc_bam.split('-')[1]
        mkdir_cmd = 'makeTagDirectory %s/%s %s/%s' % (data_dir, individual_id, data_dir, loc_bam)
        my.f_shell_cmd(mkdir_cmd)


    dir_list = my.f_grep_files_from_dir(data_dir, 'NA.*', path=True)
    my.f_print_list(dir_list)


    data_dirs=' '.join([loc_dir + '/' for loc_dir in dir_list])
    tf_dir = '%s/data/raw_data/tf/encode_peaks/processed/' % project_dir
    annotate_cmd = 'annotatePeaks.pl %s/%s hg19 -size given  -d %s -noann > %s/output.file' %(tf_dir, tf_peak[loc_tf], data_dirs, data_dir)
    my.f_shell_cmd(annotate_cmd)






#annotatePeaks.pl /homed/home/shi/expression_var/data/raw_data/tf/encode_peaks/processed/haib-gm12878-pu1.narrowPeak hg19 -size given  -d NA10851-PU1-Rep1/ NA10852-PU1-Rep1/ -noann > output.file

#RNA-seq


