

import os

home_dir = os.path.expanduser('~')
lib_dir = '%s/python/' % home_dir
import sys
sys.path.insert(0, lib_dir)
import pandas as pd
import p_mymodule as my

project_dir = '%s/expression_var/' % home_dir 

class fantom_promoter:
    file_path = None


    def __init__(self, file_path):
        self.file_path = file_path

    def extract_bed(self, col_name):
        print ''


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




print list(extract_data.loc[:,'Unnamed: 0'].str.split(':|-'))
print pd.DataFrame( list( extract_data.loc[:,'Unnamed: 0'].str.split(':|-') )).head() 
print extract_data.head()





