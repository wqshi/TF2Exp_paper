# Assign Hic fragment_id to the promoter regions.
import p_project_metadata as p_project_metadata
reload(p_project_metadata)
from p_project_metadata import *

hic_file = '%s/data/hic/high_interaction_hic.addChr.chr1.txt' % project_dir
#hic_file = '%s/data/hic/high_interaction_hic.addChr.txt.backup2' % project_dir
hic_id_pd = pd.read_csv(hic_id_file,sep='\t',header = None,names = ['chr_num', 'start', 'end', 'middle'])
chr_str = 'chr1'

print 'Shape of hic_id_pd', hic_id_pd.shape

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

log('unique_fragment_table size', unique_fragment_table.data.shape)

#2. Load the fragment-fragment interactiyon hic data, focus only one chr.
print 'Load the data'
hic_pd = pd.read_csv(hic_file, header = None, sep='\t', names = ['chr', 'cor' ,'pair_i', 'pair_j'], dtype={'pair_i':np.object, 'pair_j':np.object})
hic_pd = hic_pd[hic_pd.chr == chr_str] #select chr

#The position_i part, hic_fragment_id is the promoter id.
hic_pd.columns = ['chr', 'cor' , 'pair' ,'hic_fragment_id']

print 'Check cor score of each hic interaction ', sum(hic_pd.cor > 0 ) == hic_pd.shape[0]


unique_fragment_table.merge_feature(hic_pd, expected_cols = ['chr', 'hic_fragment_id'], check_index = False)

print unique_fragment_table.head()
print 'Check cor and pair missing at the same time:', unique_fragment_table.data.ix[:,'cor'].isnull().sum() == unique_fragment_table.data.ix[:,'pair'].isnull().sum()


print 'Length of the overlapping elements pair_i'
print len(set(unique_fragment_table.data['hic_fragment_id']).intersection(set(hic_pd['hic_fragment_id'])))

#The position_j part
hic_pd.columns = ['chr', 'cor' ,'hic_fragment_id', 'pair']
print unique_fragment_table2.head()
unique_fragment_table2.merge_feature(hic_pd, expected_cols = ['chr', 'hic_fragment_id'], check_index = False)
print 'Length of the overlapping elements pair_j'
print len(set(unique_fragment_table2.data['hic_fragment_id']).intersection(set(hic_pd['hic_fragment_id'])))
print unique_fragment_table2.data.shape
print unique_fragment_table.head()
#merge the position_i and position_j data
merge_pair_pd = pd.concat( [unique_fragment_table.data, unique_fragment_table2.data])

print merge_pair_pd.head()

unique_fragment_table.data = merge_pair_pd
unique_fragment_table.data.drop_duplicates(['chr','hic_fragment_id','pair'], inplace = True)
print unique_fragment_table.data.shape
unique_fragment_table.save_data()

print merge_pair_pd.shape
print unique_fragment_table.data.shape
print unique_fragment_table.file_path

print len( list(set(unique_fragment_table.data['pair']) - (set(unique_fragment_table2.data['hic_fragment_id']))) )

print 'Check cor and pair missing at the same time:', unique_fragment_table.data.ix[:,'cor'].isnull().sum()
print unique_fragment_table.data.ix[:,'pair'].isnull().sum()


#Step 3. Merge the hic interaction pairs back to genes.
gene_promoter_enhancer_table = region_table( "%s/assign_gene_hic_id.%s" % (output_dir, chr_str))
gene_promoter_enhancer_table.change_file_name('gene_promoter_enhancer_table' + '.' + chr_str)
gene_promoter_enhancer_table.data = gene_promoter_enhancer_table.data.drop('pair', axis = 1)
print gene_promoter_enhancer_table.data.columns
gene_promoter_enhancer_table.data.drop_duplicates(inplace = True)
gene_promoter_enhancer_table.merge_feature(unique_fragment_table.data[['chr','hic_fragment_id','pair','cor']], expected_cols = ['chr', 'hic_fragment_id'], check_index = False) 
print unique_fragment_table.data

print 'Final gene interactions after merge:', gene_promoter_enhancer_table.data.shape

print gene_promoter_enhancer_table.data.ix[:,'pair'].isnull().sum()
print gene_promoter_enhancer_table.data.ix[:,'cor'].isnull().sum()


#Step 4. Change the interaction pair (promoter-enhancer, hic_fragment_id - pair) table to one column 
print 'Columns of the gene_promoter_enhancer_talbe', gene_promoter_enhancer_table.data.columns.values
print 'Head of gene_promoter_enhancer_talbe', gene_promoter_enhancer_table.data.columns.values

gene_promoter_fragment = gene_promoter_enhancer_table.data[['chr', 'start', 'end', 'transcript_id', 'gene', 'hic_fragment_id']]
gene_promoter_fragment['cor'] = '1'
gene_promoter_fragment['type'] = 'promoter'
print gene_promoter_fragment.head()
#gene_promoter_fragment.columns = ['chr', 'start', 'end', 'transcript_id', 'gene', 'fragment']

#enhancers are at the pair column
gene_enhancer_fragment = gene_promoter_enhancer_table.data[['chr', 'start', 'end', 'transcript_id', 'gene', 'pair', 'cor' ]]
gene_enhancer_fragment['type'] = 'enhancer'
print gene_enhancer_fragment.head()

gene_enhancer_fragment = gene_enhancer_fragment[~gene_enhancer_fragment.ix[:,'cor'].isnull()]
gene_enhancer_fragment.columns = gene_promoter_fragment.columns
print gene_promoter_fragment.head()
print gene_enhancer_fragment.head()

print sum(gene_enhancer_fragment.ix[:,'cor'] == '')
print gene_enhancer_fragment.ix[:,'cor'].values

 
gene_regulatory_fragement = pd.concat([gene_promoter_fragment, gene_enhancer_fragment]).drop_duplicates()
logging.debug(gene_regulatory_fragement.head())

output_file = '%s/gene_regulatory_fragment.%s'  % (output_dir, chr_str)
gene_regulatory_fragement.to_csv(output_file, sep = '\t')

print gene_regulatory_fragement.head()


#Testing
test_fragment_id = '939687'

chr1_data = gene_promoter_enhancer_table.data[gene_promoter_enhancer_table.data.chr == 'chr1']

print sum(chr1_data.hic_fragment_id == test_fragment_id)
print ''






