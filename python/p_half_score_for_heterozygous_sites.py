#According to the genotype data, half the score if the variation is heterozygous.

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
    batch_name = 'test'
    #batch_name = '445samples_region'
    #batch_name = '462samples_sailfish'
    #batch_name = '462samples_quantile_rmNA'
batch_output_dir = f_get_batch_output_dir(batch_name)


deepsea_dir = '%s/deep_result/%s/%s/evalue/' % (batch_output_dir, mode_str ,chr_num)
vcf_dir = '%s/chr_vcf_files/%s/%s/' % (batch_output_dir, mode_str ,chr_num)



def f_half_score_het_sites(variation_file, deepsea_dir, batch_name, half_flag = True, debug = False):
    sample_id = os.path.basename(variation_file).split('.')[0]
    print sample_id
    variation_data = pd.read_csv(variation_file, sep =',', compression= 'gzip').dropna()
    if my.f_get_server_name() == 'loire' and debug == True:
        import ipdb; ipdb.set_trace()
    #print variation_data
    #variation_data.pos = variation_data.pos.astype(np.int64)
    variation_data = f_check_loc_cols(variation_data)
    
    #print variation_data.ix[1105:1110,1:6].head()
    genotype_data=variation_data[['name']]
    
    #print genotype_data.shape
    #print pd.isnull(variation_data).any(1).nonzero()[0]
    #print genotype_data[:5]
    assert len(genotype_data) == variation_data.shape[0], 'Size of genotype data not match'
    
    pattern = re.compile('(0.1|1.0)')
    het_selection = [ re.match(pattern, genotype) is not None  for genotype in genotype_data.ix[:,'name'].tolist() ]
    #import ipdb; ipdb.set_trace()
    assert len(het_selection) == variation_data.shape[0], 'Het selection is different from the variation_data'
    
    numeric_cols =my.grep_list('GM12878', variation_data.columns)
    
    het_variation = variation_data.copy()
    if half_flag == True:
        het_variation.ix[het_selection, numeric_cols] = 0.5 * variation_data.ix[het_selection, numeric_cols]
    assert all(het_variation.pos == variation_data.pos)
    assert variation_data.pos.dtype == 'int64', 'Pos data type error'
    het_variation.to_csv( "%s/%s_het/%s.diff"%(deepsea_dir, batch_name, sample_id), index = False, sep = ',', float_format='%.4e')
    

import unittest
class TestDatabaseTable(unittest.TestCase):
    def setUp(self):
        a = 0
        
    def test_half_score(self):
        variation_file_list = my.f_shell_cmd( "find %s -name '*.diff.gz'"%(deepsea_dir), quiet = True).split('\n')[0:-1]
        f_half_score_het_sites(variation_file_list[0], deepsea_dir, batch_name)
        
        half_file_list = my.f_shell_cmd( "find %s/%s_het/ -name '*.diff'"%(deepsea_dir, batch_name), quiet = True).split('\n')[0:-1]
        variation_data = pd.read_csv(variation_file_list[0], sep =',', compression= 'gzip', dtype={'pos': int}).dropna()
        het_data = pd.read_csv(half_file_list[0], sep =',').dropna()

        print len(my.grep_list('1.1',variation_data.name))

        print variation_data.pivot_table("chr",rows="name",aggfunc=len)

        
        col_name = 'GM12878|DNase|None'
        variation_data =variation_data.ix[variation_data.ix[:, col_name]!=0,:]

        print type(variation_data[[col_name]])
        print type(variation_data.ix[:, col_name])
        
        het_data = het_data.ix[het_data.ix[:, col_name]!=0,:]
        #import ipdb; ipdb.set_trace()
        
        equal_rows = abs( variation_data.ix[:, col_name] - het_data.ix[:, col_name] ) < 1e-20
        homo_rows = variation_data.name == '1|1'
        assert sum(equal_rows) == sum(homo_rows), 'Homo sites are not eaqual: not altered rows %s, home rows %s, total rows %s. This place is hard to satisfied because of approximation probelm. Exprect 0.9 overlap' %(sum(equal_rows), sum(homo_rows), len(homo_rows)  )

    def test_NA_and_float(self):
        if my.f_debug(False):
            import ipdb; ipdb.set_trace()
            
        #"The problem when the data has NA, it will convect the start to float instead the integer"
        #"This will lead to problem when use to_csv."
        deepsea_dir2 = '%s/float_test/' % deepsea_dir
        variation_file_list = my.f_shell_cmd( "find %s -name '*.diff.gz'"%(deepsea_dir2), quiet = True).split('\n')[0:-1]
        my.f_ensure_make_dir('%s/%s_het/'%(deepsea_dir2, batch_name))
        f_half_score_het_sites(variation_file_list[0], deepsea_dir2, batch_name, debug = False)
        half_file_list = my.f_shell_cmd( "find %s/%s_het/ -name '*.diff'"%(deepsea_dir2, batch_name), quiet = True).split('\n')[0:-1]
        variation_data = pd.read_csv(half_file_list[0], sep =',').dropna()
        print variation_data.columns
        assert variation_data.pos.dtype != 'float64', 'Pos data type error'
        

        
if __name__ == "__main__":
    my.f_ensure_make_dir('%s/%s_het/'%(deepsea_dir, batch_name))
    
    suite = unittest.TestLoader().loadTestsFromTestCase( TestDatabaseTable )
    unittest.TextTestRunner(verbosity=1,stream=sys.stderr).run( suite )
    
    #import ipdb; ipdb.set_trace()
    half_flag = 'homo' not in batch_name

    logging.info('Finish the test')
    print  "find %s -name '*.diff.gz'"%(deepsea_dir)
    variation_file_list = my.f_shell_cmd( "find %s -name '*.diff.gz'"%(deepsea_dir), quiet = True).split('\n')[0:-1]
    logging.info('After find')
    print 'half_flag : %s' % half_flag
    
    for loc_variation_file in variation_file_list:
        print loc_variation_file
        f_half_score_het_sites(loc_variation_file, deepsea_dir, batch_name, half_flag)

