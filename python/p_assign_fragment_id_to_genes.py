#

# Assign Hic fragment_id to the promoter regions. All the chrs.
from p_project_metadata import *
#reload(p_project_metadata)
print lib_dir
import os
import sys
sys.path.insert(0, lib_dir)
import pandas as pd
import p_mymodule as my
import pybedtools
import logging
import random
from p_region_table import region_table
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
    #batch_name = '54samples'
    batch_name = 'test'


batch_output_dir = f_get_batch_output_dir(batch_name)

#batch_name = '54samples'
#batch_name = '462samples'


gene_expression_file = f_get_batch_output_dir(batch_name) + '/rnaseq/transcript_data.bed'

hic_id_pd = pd.read_csv(hic_id_file,sep='\t',header = None,names = ['chr_num', 'start', 'end', 'middle'])
#import p_region_table as p_region_table
#reload(p_region_table)




if True:
    expression_table = region_table(gene_expression_file)
    combined_df = expression_table.data
    print 'Number of genes in the rnaseq data:', expression_table.data.shape
    expression_table.change_file_name('assign_gene_hic_id')
    print expression_table.head()
    #expression_table.data['start'] = expression_table.data['start'] - 2000
    #expression_table.save_data()
    #Clean the expression table data
    #print expression_table.data.columns

    
    feature_overlap_data = expression_table.overlap_with_feature_bed(hic_id_file, value_col=3, value_name='hic_fragment_id')
    log('Feature data head:', feature_overlap_data.head())
    print 'Size of feature overlap data:', feature_overlap_data.shape
    feature_overlap_data['hic_fragment_id'] = feature_overlap_data['hic_fragment_id'].astype(str)

    #Merge the feature based on chr, and start
    expression_table.merge_feature(feature_overlap_data)
    print 'Size of expression table', expression_table.data.shape
    print 'Number of gene with hic contact: ', expression_table.data.shape[0] - sum(expression_table.data.duplicated(['gene']))
    print 'Number of duplicates: %s' % sum(expression_table.data.duplicated())
    logging.debug('Hic_fragment_id data type:')
    print feature_overlap_data['hic_fragment_id'].dtype

    print 'Empty hic fragment ID:', sum(expression_table.data['hic_fragment_id'] == '.')

    
    #Testing random 10 samples
    i = 2
    
    combined_df.reset_index(inplace = True, drop = True)
    feature_df = feature_overlap_data
    gene_df = expression_table.data


    hic_id_df = pd.read_csv(hic_id_file, header = None, sep='\t', names = ['chr', 'start' ,'end', 'hic_fragment_id'])
    hic_id_df.index = hic_id_df.hic_fragment_id.astype(str)
    print hic_id_df.head()
    
    for i in random.sample(combined_df.index, 10):

        
        print ''
        print ''
        print '===========%s==============' % i
        #fragment_id = combined_df.ix[i, 'hic_fragment_id']
        transcript_id = combined_df.ix[i, 'gene']
    
        combined_subset = combined_df[(combined_df.gene == transcript_id) ]
        log('Combined_subset', combined_subset[['transcript_id','strand','start','end','chr']].head(10))
        #print sum(feature_df.hic_fragment_id == fragment_id)
            
        #feature_subset = feature_df[feature_df.hic_fragment_id == fragment_id]
       
  
        gene_subset = gene_df[(gene_df.gene == transcript_id)]
        log('gene_subset', gene_subset.head(10))
        print (gene_subset)
        assert all(gene_subset.start <= gene_subset.hic_fragment_id.astype(float))
        assert all(gene_subset.end >= gene_subset.hic_fragment_id.astype(float))
        
        #logging.info("Row %s, Combined %s = Gene %s * feature %s" %(i, combined_subset.shape[0], gene_subset.shape[0], feature_subset.shape[0]))
        #assert combined_subset.shape[0] == gene_subset.shape[0] * feature_subset.shape[0]




