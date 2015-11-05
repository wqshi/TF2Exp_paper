
# Assign Hic fragment_id to the promoter regions.
from p_project_metadata import home_dir, lib_dir, project_dir, hic_id_file, gene_expression_file, hic_file

import os
import sys
sys.path.insert(0, lib_dir)
import pandas as pd
import p_mymodule as my
import pybedtools
import logging



import p_region_table as p_region_table
reload(p_region_table)
from p_region_table import *

hic_file = '%s/data/hic/high_interaction_hic.addChr.txt' % project_dir
#hic_file = '%s/data/hic/high_interaction_hic.addChr.txt.backup2' % project_dir
hic_id_pd = pd.read_csv(hic_id_file,sep='\t',header = None,names = ['chr_num', 'start', 'end', 'middle'])
chr_str = 'chr1'

output_dir = '%s/data/rnaseq/' % project_dir

#1. Produce the unique list of fragment overlaped with genes.
unique_fragment_table = region_table( "%s/assign_gene_hic_id" % output_dir)
unique_fragment_table.subset_one_chr(new_file_name = 'unique',chr_str = chr_str)
print unique_fragment_table.file_path
unique_fragment_table.data.drop_duplicates(['chr', 'hic_fragment_id'], inplace= True)
print unique_fragment_table.data.shape
print unique_fragment_table.data.head(10)
unique_fragment_table2 = region_table( "%s/assign_gene_hic_id" % output_dir)
unique_fragment_table2.subset_one_chr(new_file_name = 'unique2',chr_str = chr_str)


#2. Load the fragment-fragment interaction hic data.
print 'Load the data'
hic_pd = pd.read_csv(hic_file, header = None, sep='\t', names = ['chr', 'pair_i', 'pair_j'], dtype={'pair_i':np.object, 'pair_j':np.object})
hic_pd = hic_pd[hic_pd.chr == chr_str]

#The position_i part
hic_pd.columns = ['chr', 'pair' ,'hic_fragment_id']
unique_fragment_table.merge_feature(hic_pd, expected_cols = ['chr', 'hic_fragment_id'], check_index = False)
print 'Length of the overlapping elements pair_i'
print len(set(unique_fragment_table.data['hic_fragment_id']).intersection(set(hic_pd['hic_fragment_id'])))

#The position_j part
hic_pd.columns = ['chr', 'hic_fragment_id', 'pair']
unique_fragment_table2.merge_feature(hic_pd, expected_cols = ['chr', 'hic_fragment_id'], check_index = False)
print unique_fragment_table2.file_path
print 'Length of the overlapping elements pair_j'
print len(set(unique_fragment_table2.data['hic_fragment_id']).intersection(set(hic_pd['hic_fragment_id'])))
print list(set(unique_fragment_table2.data['hic_fragment_id']).intersection(set(hic_pd['hic_fragment_id'])))[0:10]
print list(set(unique_fragment_table.data['hic_fragment_id']).intersection(set(hic_pd['hic_fragment_id'])))[0:10]
print unique_fragment_table2.data.shape


#merge the position_i and position_j data
merge_pair_pd = pd.concat( [unique_fragment_table.data, unique_fragment_table2.data])


unique_fragment_table.data = merge_pair_pd
unique_fragment_table.data.drop_duplicates(['chr','hic_fragment_id','pair'], inplace = True)
print unique_fragment_table.data.shape
unique_fragment_table.save_data()

print merge_pair_pd.shape
print unique_fragment_table.data.shape
print unique_fragment_table.file_path

print len( list(set(unique_fragment_table.data['pair']) - (set(unique_fragment_table2.data['hic_fragment_id']))) )



#Step 3. Merge the hic interaction pairs back to genes.
gene_promoter_enhancer_table = region_table( "%s/assign_gene_hic_id.%s" % (output_dir, chr_str))
gene_promoter_enhancer_table.change_file_name('gene_promoter_enhancer_table' + '.' + chr_str)
gene_promoter_enhancer_table.merge_feature(unique_fragment_table.data[['chr','hic_fragment_id','pair']], expected_cols = ['chr', 'hic_fragment_id'], check_index = False) 
print unique_fragment_table.data
print 'Final gene interactions after merge:', gene_promoter_enhancer_table.data.shape




#Testing
test_fragment_id = '939687'

chr1_data = gene_promoter_enhancer_table.data[gene_promoter_enhancer_table.data.chr == 'chr1']

print sum(chr1_data.hic_fragment_id == test_fragment_id)
print ''













