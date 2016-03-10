

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
    #batch_name = '462samples_sailfish'
    #batch_name = '462samples_quantile_rmNA'
batch_output_dir = f_get_batch_output_dir(batch_name)

tf_dir = '%s/data/raw_data/tf/encode_peaks/' % project_dir
deepsea_dir = '%s/deep_result/%s/%s/evalue/' % (batch_output_dir, mode_str ,chr_num)
vcf_dir = '%s/chr_vcf_files/%s/%s/' % (batch_output_dir, mode_str ,chr_num)



def f_half_score_het_sites(variation_file, deepsea_dir):
    sample_id = os.path.basename(variation_file).split('.')[0]
    variation_data = pd.read_csv(variation_file, sep =',', compression= 'gzip').dropna()
    print sample_id
    #print variation_data.ix[1105:1110,1:6].head()
    genotype_data=variation_data[['name']]

    #print genotype_data.shape
    #print pd.isnull(variation_data).any(1).nonzero()[0]
    #print genotype_data[:5]
    assert len(genotype_data) == variation_data.shape[0], 'Size of genotype data not match'

    
    pattern = re.compile('(0.1|1.0)')
    het_selection = [ re.match(pattern, genotype) is not None  for genotype in genotype_data.ix[:,'name'].tolist() ]

    assert len(het_selection) == variation_data.shape[0], 'Het selection is different from the variation_data'

    numeric_cols =my.grep_list('GM12878', variation_data.columns)
    
    het_variation = variation_data.copy()
    het_variation.ix[het_selection, numeric_cols] = 0.5 * variation_data.ix[het_selection, numeric_cols]

    het_variation.to_csv( "%s/het/%s.diff"%(deepsea_dir, sample_id), index = False, sep = ',')



###Part 2. Correct the bias of GM12878





import unittest
class TestDatabaseTable(unittest.TestCase):
    def setUp(self):
        a = 0
        
    def test_half_score(self):
        variation_file_list = my.f_shell_cmd( "find %s -name '*.diff.gz'"%(deepsea_dir), quiet = True).split('\n')[0:-1]

        
        f_half_score_het_sites(variation_file_list[0], deepsea_dir)
        
        half_file_list = my.f_shell_cmd( "find %s/het/ -name '*.diff'"%(deepsea_dir), quiet = True).split('\n')[0:-1]
        variation_data = pd.read_csv(variation_file_list[0], sep =',', compression= 'gzip').dropna()
        het_data = pd.read_csv(half_file_list[0], sep =',').dropna()

        print len(my.grep_list('1.1',variation_data.name))

        print variation_data.pivot_table("chr",rows="name",aggfunc=len)

        
        col_name = 'GM12878|DNase|None'
        variation_data =variation_data.ix[variation_data.ix[:, col_name]!=0,:]

        print type(variation_data[[col_name]])
        print type(variation_data.ix[:, col_name])
        
        het_data = het_data.ix[het_data.ix[:, col_name]!=0,:]
        print variation_data.columns
        print variation_data.shape
        
        
        
        equal_rows = abs( variation_data.ix[:, col_name] - het_data.ix[:, col_name] ) < 0.0000000000000001
        homo_rows = variation_data.name == '1|1'
        assert sum(equal_rows) == sum(homo_rows), 'Homo sites are not eaqual by %s out of %s' %(sum(equal_rows) - sum(homo_rows), len(homo_rows)  )








    
if __name__ == "__main__":
    my.f_ensure_make_dir('%s/het/'%deepsea_dir)
    
    suite = unittest.TestLoader().loadTestsFromTestCase( TestDatabaseTable )
    unittest.TextTestRunner(verbosity=1,stream=sys.stderr).run( suite )

    
    variation_file_list = my.f_shell_cmd( "find %s -name '*.diff.gz'"%(deepsea_dir), quiet = True).split('\n')[0:-1]

    for loc_variation_file in variation_file_list:
        f_half_score_het_sites(loc_variation_file, deepsea_dir)

