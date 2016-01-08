
#extract the deepsea predictions of one tf.
#Deepsea prediction is sample based. The output of this script is TF based.

from p_project_metadata import *

import re
import tempfile
#reload(p_region_table)

import argparse

parser = argparse.ArgumentParser(description='Extract the deepsea predictions of one tf.Deepsea prediction is sample based. The output of this script is TF based.')


print '==========',__doc__

if __doc__ is None:
    parser.add_argument("--chr",     help="The chr wanted to compute", default="chr22")
    parser.add_argument("--mode", help="Output type: test(only first 1000 variants); all(all the variants)", default="all")
    parser.add_argument('--batch_name', help = "462samples or 54samples", default = '462samples' )
    args = parser.parse_args()
    chr_num = args.chr
    mode_str = args.mode
    batch_name = args.batch_name
else:
    chr_num = 'chr22'
    mode_str = 'all'
    batch_name = '54samples'
batch_output_dir = f_get_batch_output_dir(batch_name)

tf_dir = '%s/data/raw_data/tf/encode_peaks/' % project_dir
deepsea_dir = '%s/deep_result/%s/%s/' % (batch_output_dir, mode_str ,chr_num)

#deepsea_data = pd.read_csv('%s/infile.head.txt'%deepsea_dir, sep = ',')

peak_list_raw = my.f_shell_cmd( "find %s -name '*gm12878-*.narrowPeak'"%(tf_dir), quiet = True).split('\n')

black_list = my.grep_list(".*(--|Rep[1-9]|-myc|-pu1)", peak_list_raw)
peak_list = list(set(peak_list_raw) - set(['']) - set(black_list))




tf_list = [[re.split('[.]|-',os.path.basename(peak_file))[2], os.path.getsize(peak_file)] for peak_file in peak_list ]

peak_file_df = pd.DataFrame(data = tf_list, columns = ['tf', 'file_size'])

peak_file_df['file_path'] = peak_list


##Keep the large peak file.
peak_file_df_rmdup = peak_file_df.sort(columns = 'file_size', ascending = False ).drop_duplicates(['tf'])
print peak_file_df_rmdup.shape
print peak_file_df.ix[1:3,:]
print peak_file_df.ix[1:3,:].sort(columns = 'file_size', ascending = False ).drop_duplicates(['tf'])
print peak_file_df_rmdup.head()





############################

loc_tf_list = ['ctcf','ebf1']
variation_file_list = my.f_shell_cmd( "find %s -name '*.diff'"%(deepsea_dir), quiet = True).split('\n')[0:-1]

variation_data = pd.read_csv(variation_file_list[0], sep =',')
loc_tf_list = list(set( [ feature.split('|')[1].lower() for feature in my.grep_list('GM12878',variation_data.columns.values)]))
my.f_print_list(loc_tf_list)
logging.info('Number of features %s' % len(loc_tf_list) )

tf_variation_dir = '%s/output/tf_variation/%s/%s' % (batch_output_dir, mode_str, chr_num)
my.f_ensure_make_dir(tf_variation_dir)
target_cell = 'gm12878'





for loc_tf in loc_tf_list:
#for loc_tf in ['c-fos']:
    try:
        tf_peak_file = peak_file_df_rmdup.ix[peak_file_df_rmdup.tf == loc_tf, 'file_path'].tolist()[0]
    except:
        logging.error("Don't find the peak file for %s " % loc_tf)
        continue
    print tf_peak_file
    binding_regions = pd.io.parsers.read_csv(tf_peak_file, sep="\t", header=None).ix[:,0:2]
    
    binding_regions.columns = ['chr', 'start', 'end']
    binding_regions = binding_regions[binding_regions.chr == chr_num]
    print binding_regions.head()
    
    tf_regions_table = region_table('%s/%s_matrix.txt'%(tf_variation_dir, loc_tf), binding_regions)

    print tf_regions_table.file_path
    
    for variation_file in variation_file_list:
        sample_id = os.path.basename(variation_file).replace('.diff', '')
        variation_data = pd.read_csv(variation_file, sep =',')
        variation_data['start'] = variation_data['pos']
        variation_data['end'] =  variation_data['pos'] + 1
        #print variation_table.head()
        
        #extract_variation_data from variation_table
        tf_variation_cols = my.grep_list('%s.*%s'%(target_cell, loc_tf), variation_data.columns)

        if (len(tf_variation_cols) == 0):
            logging.warning('Find 0 matched prediction for %s'%loc_tf)
        elif (len(tf_variation_cols) == 0) > 1 :
            logging.warning('Find more than 1 matched prediction for %s'%loc_tf)

        #tmp_bed_file = tempfile.NamedTemporaryFile(dir = tf_variation_dir, prefix = 'tmp.bed' , delete=False, bufsize = 10000000)

        tmp_bed_file = tf_variation_dir + '/' + my.f_generate_tmp_file_name('bed')

        tf_selected_column = tf_variation_cols[0]
        
        variation_data.ix[:,['chr', 'start', 'end', tf_selected_column ]].to_csv(tmp_bed_file , header = False, index = False ,sep = '\t')
        #print variation_data.columns
        #print variation_data.shape
        
        #Intersect with the tf_binding regions.
        tf_variation_data = tf_regions_table.overlap_with_feature_bed(tmp_bed_file, 3, value_name=sample_id)
        #print tf_regions_table.data.shape
        
        #print tf_variation_data.head()
        
        #Aggregete the impact for the same regions
        tf_variation_data[sample_id] = tf_variation_data[sample_id].astype(float)        
        agg_variation_data = tf_variation_data.groupby(['chr', 'start'], as_index=False).sum()

        #merge aggregated impact to the TF binding regions.
        tf_regions_table.merge_feature(agg_variation_data)

        os.remove(tmp_bed_file)


    os.remove(tf_regions_table.loc_file)

tf_file_list = my.f_shell_cmd( "find %s -name '*_matrix.txt'"%(tf_variation_dir), quiet = True).split('\n')[0:-1]
tf_name_list = [ os.path.basename(tf_file).replace('_matrix.txt', '').upper() for tf_file in tf_file_list]
my.f_print_list(tf_file_list)
my.f_print_list(tf_name_list)

tf_summary_list  = pd.DataFrame([tf_name_list, tf_file_list]).T
tf_summary_list.columns = ['feature', 'path']
tf_summary_list.to_csv('%s/data_path.csv'%tf_variation_dir, index = False, sep = ' ')

print tf_summary_list.shape

