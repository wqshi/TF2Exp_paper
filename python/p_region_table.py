from p_project_metadata import home_dir, lib_dir
import os
import sys
sys.path.insert(0, lib_dir)
import pandas as pd
import p_mymodule as my
import pybedtools
import logging
import p_pd as p_pd
reload(p_pd)
from p_pd import data_table
from lockfile import LockFile
from lockfile import LockTimeout


class region_table(data_table):
    file_path = None
    data=None
    chr_str = ''
    loc_file = None
    
    def __init__(self, file_path, data=None):
        data_table.__init__(self, file_path=file_path, data=data)
        #self.data.set_index(keys=["chr","start"], inplace=True, drop=False)
        self.data = self.data.sort(['chr','start'])
        self.data.sort_index(axis=1, ascending=False, inplace=True)
        if 'hic_fragment_id' in self.data.columns:
            self.data.hic_fragment_id = self.data.hic_fragment_id.astype(str)
        #self.data.reindex_axis(sorted(self.data.columns, reverse = True),axis = 1, copy = False)
        #self.save_data()
    def change_file_name(self, new, chr_str=''):
        #import ipdb; ipdb.set_trace()
        if chr_str != '':
            chr_str = '.' + chr_str
        self.file_path = "%s/%s%s" %( os.path.dirname(self.file_path), new, chr_str)
        self.save_data()

    def subset_one_chr(self, new_file_name = None, chr_str=''):
        if chr_str != '':
            self.data = self.data[self.data.chr == chr_str]
            logging.info('After subset to %s, %s entries left' % (chr_str, self.data.shape[0]))
            
        if new_file_name is None:
            new_file_name = os.path.basename(self.file_path)
        self.change_file_name(new_file_name, chr_str)
    def save_data(self):
        self.data.drop_duplicates(inplace = True)
        self.data.to_csv(self.file_path, index=None, sep="\t", na_rep=".", float_format='%.4e')
        
    def extract_bed(self):
        #print self.data.head()
        #import ipdb; ipdb.set_trace()
        if self.loc_file is not None:
            return self.loc_file
        logging.info('Generate new bed')
        data_dir = os.path.dirname(self.file_path)
        bed_file = data_dir + '/' + my.f_generate_tmp_file_name("loc.bed")
        bed_data=self.data[["chr","start",'end']]
        bed_data["name"]=bed_data["chr"]+"-"+bed_data["start"].map(str)
        bed_data.drop_duplicates().to_csv(bed_file, header=False, index=False, sep="\t")
        self.loc_file = bed_file
        return bed_file

    def head(self):
        sample_cols = my.grep_list('(NA|HG)[0-9]+', self.data.columns)
        show_cols = list(set(self.data.columns) - set(sample_cols))
        self.data = self.data.drop(sample_cols, axis = 1)
        print self.data.ix[:,show_cols].head()
        
    def overlap_with_feature_bed(self, feature_file, value_col, value_name, feature_extend_distanace = 0, debug= False):
        #Extend the feature_file by feature_extend_distance, and overlap with loc_file, and get the overlap data
        #loc_file: the bed extracted from the database fiile, point positions
        #value_col: the column of the desired value in the feature file, 0 based
        #value_name: the name of the value_col
        #feature_extend_distance: this is only for the point features. for region features, set it to 0.
        #return: [loc_chr, loc_start, value_col], loc_start is 1 based same as in database file
        #import ipdb; ipdb.set_trace()
        import socket
        server_name = socket.gethostname()
        if debug == True and 'loire' in server_name :
            import ipdb; ipdb.set_trace()

        loc_file = self.extract_bed()
        bed_data = pybedtools.BedTool(loc_file)
            
        feature_extend_bed = feature_file
        feature_regions=bed_data.intersect(feature_extend_bed,wo=True)
        
        feature_regions_pd=my.f_bed_to_pd(feature_regions).ix[:,[0,1, 7]] # the 7th position is for value_col
        feature_regions_pd.columns=["chr","start",value_name]
        #logging.debug(feature_regions_pd[['start', value_name]].values)
        #logging.debug("=====Missing Pvalue======") #Mostly because I extend 100bp
        feature_regions_pd['start'] = feature_regions_pd['start'].astype(float)
        
        #os.remove(loc_file)
        return feature_regions_pd


    def overlap_with_feature_bed_extend(self, feature_file, value_col, value_name, feature_extend_distance = 0, debug= False):
        #Extend the feature_file by feature_extend_distance, and overlap with loc_file, and get the overlap data
        #loc_file: the bed extracted from the database fiile, point positions
        #value_col: the ith column of the desired value in the feature file, 1 based
        #value_name: the name of the value_col 
        #feature_extend_distance: this is only for the point features. for region features, set it to 0.
        #return: [loc_chr, loc_start, value_col], loc_start is 1 based same as in database file
        #import ipdb; ipdb.set_trace()
        import socket
        server_name = socket.gethostname()
        if debug == True and 'loire' in server_name :
            import ipdb; ipdb.set_trace()

        loc_file = self.extract_bed()
        bed_data = pybedtools.BedTool(loc_file)
        bed_data_extended = bed_data.slop(genome= 'hg19', b= feature_extend_distance)
        
        feature_extend_bed = feature_file
        feature_regions=bed_data_extended.intersect(feature_extend_bed,wo=True)
        
        feature_regions_pd=my.f_bed_to_pd(feature_regions).ix[:,[0,1,2,3, value_col + 3 ]]
        
        feature_regions_pd.columns=["chr","start",'end', 'name' ,value_name]
        #logging.debug(feature_regions_pd[['start', value_name]].values)
        #logging.debug("=====Missing Pvalue======") #Mostly because I extend 100bp
        feature_regions_pd['start'] = feature_regions_pd['start'].astype(float) + feature_extend_distance
        feature_regions_pd['end'] = feature_regions_pd['end'].astype(float) - feature_extend_distance
        del bed_data
        del bed_data_extended
        #os.remove(loc_file)
        return feature_regions_pd

    def delete_file(self):
        os.remove(self.file_path)
        
    def merge_feature(self, feature_data,expected_cols=["chr", "start"], check_index = True ,debug = False):

        if debug == True:
            import ipdb; ipdb.set_trace()
        #logging.debug(self.show_size(return_flag =True))
        assert set(expected_cols) < set(feature_data.columns), "Unexpected col names in feature data"
        feature_data.set_index(keys=expected_cols, inplace=True, drop=False)
        self.data.set_index(keys=expected_cols, inplace=True, drop=False)
        #print self.data.index[1:10]
        
        if check_index ==True and (set(self.data.index) >= set(feature_data.index)) == False:
            #print self.data.index
            #print feature_data.index
            #selection = 1:10
            #print set(feature_data.ix[selection,"chr"]+feature_data.ix[selection, "start"].map(str)) - set(self.data.ix[selection,"chr"]+self.data.ix[selection,"start"].map(str))
            assert set(self.data.index) >= set(feature_data.index), "Unexpected index"
  
        new_cols = list( set(feature_data.columns) - set(self.data.columns) )
        update_cols = list(set(feature_data.columns).intersection(self.data.columns) - set(expected_cols) - set(['ref']) )
        
        if  update_cols !=[]:
            duplicated_rows = my.f_duplicated_index(feature_data, expected_cols)

            if any(duplicated_rows) > 0:
                logging.info("=======duplicated data===========")
                print feature_data.ix[duplicated_rows, :]

            
            update_data = feature_data.drop_duplicates(expected_cols).ix[:,update_cols]
            #update_data = feature_data.drop_duplicates().ix[:,update_cols]
            #print update_data.ix[0:10,:]
            self.data.update(update_data)
            print 'Updated data in when merge'
            print self.data.ix[0:10,update_cols]
            #print self.data.index[1:10]
        if new_cols !=[]:
            self.data = pd.merge(self.data, feature_data.ix[:,new_cols + expected_cols], on=expected_cols, how="left")

        #self.data = tmp_data
        self.data.set_index(keys=["chr","start"], inplace=True, drop=False)
        self.save_data()
        self.show_size(return_flag = True)

