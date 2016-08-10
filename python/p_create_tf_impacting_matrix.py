
#extract the deepsea predictions of one tf.
#Deepsea prediction is sample based. The output of this script is TF based.

from p_project_metadata import *

import re
import tempfile
#reload(p_region_table)

from joblib import Parallel, delayed  
import multiprocessing
num_cores = multiprocessing.cpu_count()-2
print num_cores
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
    batch_name = 'test'
batch_output_dir = f_get_batch_output_dir(batch_name)


deepsea_dir = '%s/deep_result/%s/%s/evalue/het/' % (batch_output_dir, mode_str ,chr_num)

#deepsea_data = pd.read_csv('%s/infile.head.txt'%deepsea_dir, sep = ',')
peak_list = f_get_tf_peak_list(project_dir)

tf_list = [[re.split('[.]|-',os.path.basename(peak_file))[2], os.path.getsize(peak_file)] for peak_file in peak_list ]
peak_file_df = pd.DataFrame(data = tf_list, columns = ['tf', 'file_size'])
peak_file_df['file_path'] = peak_list

##Keep the large peak file.
peak_file_df_rmdup = peak_file_df.sort(columns = 'file_size', ascending = False ).drop_duplicates(['tf'])
print peak_file_df_rmdup.shape
print peak_file_df.ix[1:3,:]
print peak_file_df.ix[1:3,:].sort(columns = 'file_size', ascending = False ).drop_duplicates(['tf'])
print peak_file_df_rmdup.head()

#the variation distacne to tf. Not all the variants, only the variants near the Tf binding regions.
variation_dis = pd.read_csv('%s/data/raw_data/wgs/1kg/variant_tf_distance/variation_dis_to_peakMax.%s' % (project_dir, chr_num), sep = '\t',na_values = '.')
variation_dis.drop_duplicates('start', inplace = True)
variation_dis['min_dis']=variation_dis.ix[:,4:].apply(np.nanmin, axis = 1)

print len(variation_dis.columns.values[4:-1])

print set(peak_file_df_rmdup.tf) - set(variation_dis.columns.values[4:-1])


print variation_dis.shape


############################

loc_tf_list = ['ctcf','ebf1']
variation_file_list = my.f_shell_cmd( "find %s -name '*.diff'"%(deepsea_dir), quiet = True).split('\n')[0:-1]

variation_data = pd.read_csv(variation_file_list[0], sep =',')
loc_tf_list = list(set( [ feature.split('|')[1].lower() for feature in my.grep_list('GM12878',variation_data.columns.values)]))
loc_tf_list.sort()
my.f_print_list(loc_tf_list)
if batch_name == 'test':
    loc_tf_list = ['pol3','c-fos'] + loc_tf_list[1:3]
    
logging.info('Number of features %s' % len(loc_tf_list) )

tf_variation_dir = '%s/output/tf_variation/%s/%s' % (batch_output_dir, mode_str, chr_num)
my.f_ensure_make_dir(tf_variation_dir)
target_cell = 'gm12878'

def f_convert_deepseaTF_to_peakTF(loc_tf_name):
    return loc_tf_name.replace('-','').replace('.','')


#Manually 
tf_name_pd = pd.DataFrame(data = [loc_tf_list, loc_tf_list]).T
tf_name_pd.columns = ['tf', 'rename']
tf_name_pd.to_csv('%s/data/raw_data/tf_name_match.txt'%project_dir, index = False, sep='\t' )
print tf_name_pd.head()

def f_subet_variation_according_to_tf_distance(variation_data_raw, tf_name, variation_dis):
    #print ''
    #Subset or selcted the variations to which are closest to the peak max of targeted TF.
    variation_dis.index  = variation_dis.start.astype(str)
    variation_data_raw.index = variation_data_raw.start.astype(str)
    intersect_pos = list(set(variation_dis.index).intersection(set(variation_data_raw.index)))
    subset_dis=variation_dis.ix[intersect_pos,:]
    assigned_locations = np.logical_and( subset_dis[tf_name] == subset_dis.min_dis, ~subset_dis.min_dis.isnull())
    variation_assigned_to_tf = variation_data_raw.ix[ subset_dis.index[assigned_locations],:]
    return variation_assigned_to_tf


def parse_one_tf(variation_file_list, peak_file_df_rmdup, target_cell, loc_tf, tf_variation_dir):
    #import ipdb; ipdb.set_trace()
    loc_peak_tf = f_convert_deepseaTF_to_peakTF(loc_tf)
    try:
        tf_peak_file = peak_file_df_rmdup.ix[peak_file_df_rmdup.tf == loc_peak_tf, 'file_path'].tolist()[0]
    except:
        logging.error("Don't find the peak file for %s " % loc_tf)
        return 0
    #import ipdb; ipdb.set_trace()
    if loc_peak_tf not in variation_dis.columns.values:
        logging.info('Missing %s in variation dis', loc_peak_tf)
        return 0

    
    print tf_peak_file
    binding_regions = pd.io.parsers.read_csv(tf_peak_file, sep="\t", header=None).ix[:,0:2]
    
    binding_regions.columns = ['chr', 'start', 'end']
    binding_regions = binding_regions[binding_regions.chr == chr_num]
    print binding_regions.head()
    
    tf_regions_table = region_table('%s/%s_matrix.txt'%(tf_variation_dir, loc_tf), binding_regions)

    print tf_regions_table.file_path
    
    tf_regions_table.extract_bed()

    
    #import ipdb; ipdb.set_trace()
    
    for variation_file in variation_file_list:
    
        sample_id = os.path.basename(variation_file).split('.')[0]
        variation_data_raw = pd.read_csv(variation_file, sep =',')
        variation_data_raw['start'] = variation_data_raw['pos'].astype(int) -1
        variation_data_raw['end'] =  variation_data_raw['pos'].astype(int)
        if 'nearest' in batch_name:
            variation_data = f_subet_variation_according_to_tf_distance(variation_data_raw, loc_peak_tf, variation_dis)
            logging.info()
        else:
            variation_data = variation_data_raw

        #print variation_data.ix[:,0:10].head()

        #extract_variation_data from variation_table
        tf_variation_cols = my.grep_list('%s.*%s'%(target_cell, loc_tf), variation_data.columns)

        if (len(tf_variation_cols) == 0):
            logging.warning('Find 0 matched prediction for %s'%loc_tf)
        elif (len(tf_variation_cols) >1):
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

    if tf_regions_table.loc_file is not None:
        os.remove(tf_regions_table.loc_file)


for loc_tf in loc_tf_list[2:]:
    print loc_tf
    parse_one_tf(variation_file_list, peak_file_df_rmdup, target_cell, loc_tf, tf_variation_dir)
    break 

Parallel(n_jobs=4)(delayed(parse_one_tf)(variation_file_list, peak_file_df_rmdup, target_cell, loc_tf, tf_variation_dir) for loc_tf in loc_tf_list)  




#Delete the half data dir.
import shutil
shutil.rmtree(deepsea_dir)



tf_file_list = my.f_shell_cmd( "find %s -name '*_matrix.txt'"%(tf_variation_dir), quiet = True).split('\n')[0:-1]
tf_name_list = [ os.path.basename(tf_file).replace('_matrix.txt', '').upper() for tf_file in tf_file_list]
my.f_print_list(tf_file_list)
my.f_print_list(tf_name_list)

tf_summary_list  = pd.DataFrame([tf_name_list, tf_file_list]).T
tf_summary_list.columns = ['feature', 'path']
tf_summary_list.to_csv('%s/data_path.csv'%tf_variation_dir, index = False, sep = ' ')

print tf_summary_list.shape




