#! python 2.7

import os

home_dir = os.path.expanduser('~')
lib_dir = '%s/python/' % home_dir
import sys
sys.path.insert(0, lib_dir)
import pandas as pd
import p_mymodule as my

project_dir = '%s/expression_var/' % home_dir
#Parse the enhancer file

loc_cell = 'gm12878'
enhancer_path = '%s/data/fantom5/hg19_permissive_enhancers_expression_rle_tpm.csv' % project_dir
enhancer_pd = pd.read_csv(enhancer_path, header = 0, sep=',', )
print enhancer_pd.shape
target_cel_columns = my.grep_list( '.*gm12878', enhancer_pd.columns.tolist())
print ['Unnamed: 0'] + target_cel_columns
extract_data = enhancer_pd.loc[:, ['Unnamed: 0'] + target_cel_columns]
extract_data
print extract_data.head()



coord_data = pd.DataFrame( list( extract_data.loc[:,'Unnamed: 0'].str.split(':|-') ) )

extract_data['chr'] = coord_data[0]
extract_data['start'] = coord_data[1]
extract_data['end'] = coord_data[2]

extract_data.columns = ['name', 'rep1', 'rep2', 'rep3', 'chr', 'start', 'end']
bed_data = extract_data.loc[:,['chr', 'start', 'end', 'rep1', 'rep2', 'rep3']]
print bed_data.head()


bed_data.to_csv('%s/data/fantom5/enhancer_%s.bed' % (project_dir, loc_cell), sep = '\t', header = False, index = False, na_rep = '.')




####promoter parse######

promoter_file = '%s/data/fantom5/hg19.cage_peak_tpm_ann_decoded.osc.txt.gz.extract.tsv' % project_dir


promoter_pd = pd.read_csv(promoter_file, header = 0, sep='\t')

promoter_pd.columns = ['name', 'short_description', 'description','transcript','entrezgene_id', 'hgnc_id', 'rep1', 'rep2', 'rep3' ]


print promoter_pd.head()

coord_data = pd.DataFrame( list( promoter_pd.loc[:,'name'].str.split(':|-|[.][.]|,') ))
promoter_pd['chr'] = coord_data[0]
promoter_pd['start'] = coord_data[1]
promoter_pd['end'] = coord_data[2]

bed_data = promoter_pd.loc[:,['chr', 'start', 'end', 'transcript' ,'rep1', 'rep2', 'rep3']]
print bed_data.head()


bed_data.to_csv('%s/data/fantom5/promoter_%s.bed' % (project_dir, loc_cell), sep = '\t', header = False, index = False,na_rep = '.')










