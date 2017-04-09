#Step1 Add the coordinate to the gene expression file.

chr_str=$2
mode='all'
#batch_name='54samples_peer'
batch_name=$1


Rscript3 ../R/s_tf_expression.R --batch_name $batch_name --sep space

echo ''
echo ''
echo 'Step 2. Assign hic_id to the genes'
python2.7 p_assign_fragment_id_to_genes.py --chr $chr_str --mode all --batch_name $batch_name



echo ''
echo ''
echo 'Step 3. Assign interacted hic_id to genes promoters.'
python2.7 p_assign_fragment_id_to_regions.py  --chr $chr_str --mode all --batch_name $batch_name


#Step 4. Create TF variation data
#python2.7 p_half_score_for_heterozygous_sites.py --chr $chr_str --mode all --batch_name $batch_name
#python2.7 p_create_tf_impacting_matrix.py           --chr $chr_str --mode all --batch_name $batch_name


echo ''
echo ''
echo 'Step 5. Add the histone and TF data.'
python2.7 p_create_gene_regulaory_whole_profile.py --chr $chr_str --mode all --batch_name $batch_name


#Step 6. Split each gene into a seperate file
Rscript3 ../R/r_split_genes.R --batch_name $batch_name --chr $chr_str

#Step 7. Regression model in R.
sh ../R/s_qsub_jobs.sh $batch_name > ./data/$batch_name/output/gene_regression.log 2>&1


