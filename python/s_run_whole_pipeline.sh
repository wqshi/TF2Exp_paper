#Step1 Add the coordinate to the gene expression file.

cd ~/expression_var/python/

chr_str=$2
mode='all'
#batch_name='54samples_peer'
batch_name=$1


###Test cases#####
##sh s_run_whole_pipeline.sh test chr22


#Rscript3 ../R/r_split_genes.R --batch_name $batch_name --chr $chr_str
#exit 0
##########

#Here can be conflicted because multiple chrs try to write the same thing.
#Rscript3 ../R/s_tf_expression.R --batch_name $batch_name --sep space


echo ''
echo ''
echo 'Step 2. Assign hic_id to the genes'
python2.7 p_assign_fragment_id_to_genes.py --chr $chr_str --mode all --batch_name $batch_name



echo ''
echo ''
echo 'Step 3. Assign interacted hic_id to genes promoters.'
python2.7 p_assign_fragment_id_to_regions.py  --chr $chr_str --mode all --batch_name $batch_name


#Step 4. Create TF variation data
echo ''
echo ''
#This step is only needed when update the deepsea results or change the size of sampels
#For instance, if we want to test a new batch e.g 462samples_new, we can use the data from 462samples_sailfish.
echo 'Step 4.1 Half score at heterozygous sites'
python2.7 p_half_score_for_heterozygous_sites.py --chr $chr_str --mode all --batch_name $batch_name

echo 'Step 4.2 [Only for the nearest model] Assign the variants to the nearest TF binding peaks'
#python2.7 p_find_closest_peak_for_variation.py --chr $chr_str --batch_name $batch_name

echo 'Step 4.3 Create TF variation data'
python2.7 p_create_tf_impacting_matrix.py           --chr $chr_str --mode all --batch_name $batch_name

#exit 0
echo ''
echo ''
echo 'Step 5. Add the histone and TF data.'
python2.7 p_create_gene_regulaory_whole_profile.py --chr $chr_str --mode all --batch_name $batch_name


#Step 6. Split each gene into a seperate file
echo 'Step 6. Split genes'
Rscript3 ../R/r_split_genes.R --batch_name $batch_name --chr $chr_str


#Step 7. Regression model in R. Run every thing in loire
#sh ../R/s_qsub_jobs.sh $batch_name > ./data/$batch_name/output/gene_regression.log 2>&1


#Step 8. Run the regression model in clustdell
#sh ../R/s_start_cluster_gene_job.sh $batch_name enet $chr_str


