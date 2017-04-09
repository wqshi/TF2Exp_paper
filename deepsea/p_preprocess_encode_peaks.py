import os
home_dir = os.path.expanduser('~')
lib_dir = '%s/python/' % home_dir
import sys
sys.path.insert(0, lib_dir)
sys.path.insert(0, '%s/expression_var/python/' % home_dir)
import pandas as pd
import p_mymodule as my
from p_project_metadata import *
from p_generate_peak_fastq import chipseq_region

peak_file_df_rmdup = f_get_peak_file_df_rmdup(project_dir)

print peak_file_df_rmdup.head()

processed_dir = '%s/data/raw_data/tf/encode_peaks/processed/' % project_dir
my.f_ensure_make_dir(processed_dir)



for loc_tf in peak_file_df_rmdup.tf:
    #loc_tf = 'pol2'
    print 'Process %s' % loc_tf
    peak_file = peak_file_df_rmdup.ix[loc_tf, 'file_path']
    tf_region = chipseq_region(file_path = peak_file)
    tf_region.merge_overlapped_peaks()
    tf_region.split_peaks_with_multiple_peakMax(debug = False)
    print tf_region.binding_df.head()
    #import ipdb; ipdb.set_trace()
    #print tf_region.binding_df.ix[ tf_region.binding_df.start== 43044464,:]
    tf_region.bed_trim_binding_regions()
    tf_region.save_bed('%s/%s' %(processed_dir, os.path.basename(peak_file)))
