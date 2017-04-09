####################Using 462sample_log as base to fit the 462_peer model ####################################
batch_name='462samples_genebody'
test_flag='FALSE' #binary
gene='gene'
add_miRNA='FALSE'
add_TF_exp='FALSE'
add_permutation='FALSE'
add_TF_exp_only='FALSE'
add_predict_tf='FALSE'
add_TF_exp_only='FALSE'
add_YRI='TRUE'
chr_str='chr22'
population='all'
TF_exp_type='fakeTF'
model='glmnet'
add_gm12878='TRUE'
new_batch='462samples_genebody'


#chr_array=({22..1})
#chr_array=(22 20 10 2)
chr_array=(22)
for chr_num in ${chr_array[@]}
do
    chr_str='chr'$chr_num
    new_batch_random='all'
    #sh ./s_start_cluster_gene_job.sh $batch_name $test_flag $model $chr_str $gene $add_miRNA $add_TF_exp $add_permutation $add_TF_exp_only $add_predict_tf $add_YRI $population $TF_exp_type $add_gm12878 $new_batch $new_batch_random

    new_batch_random='TF'
    sh ./s_start_cluster_gene_job.sh $batch_name $test_flag $model $chr_str $gene $add_miRNA $add_TF_exp $add_permutation $add_TF_exp_only $add_predict_tf $add_YRI $population $TF_exp_type $add_gm12878 $new_batch $new_batch_random

    new_batch_random='SNP'
    #sh ./s_start_cluster_gene_job.sh $batch_name $test_flag $model $chr_str $gene $add_miRNA $add_TF_exp $add_permutation $add_TF_exp_only $add_predict_tf $add_YRI $population $TF_exp_type $add_gm12878 $new_batch $new_batch_random

    new_batch_random='random'
    #sh ./s_start_cluster_gene_job.sh $batch_name $test_flag $model $chr_str $gene $add_miRNA $add_TF_exp $add_permutation $add_TF_exp_only $add_predict_tf $add_YRI $population $TF_exp_type $add_gm12878 $new_batch $new_batch_random
done

exit 0
