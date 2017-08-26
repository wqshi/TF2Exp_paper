batch_name=$1
test_flag=$2
model=$3
chr_str=$4
gene=$5
add_miRNA=$6
add_TF_exp=$7
add_permutation=$8
add_TF_exp_only=$9
add_predict_tf=${10}
add_YRI=${11}
population=${12}
TF_exp_type=${13}
add_gm12878=${14}
new_batch=${15}
new_batch_random=${16}
project_dir=$HOME'/expression_var/R/'
if [ $2 != 'classif' ]; then
    echo 'Regression'
    Rscript3 $project_dir/s_regression_for_one_gene_OOB.R --batch_name $batch_name --add_histone FALSE \
                   --test FALSE --model $model --chr_str $chr_str --gene $gene --add_miRNA $add_miRNA  \
                   --add_TF_exp $add_TF_exp --add_penalty $add_permutation --add_TF_exp_only $add_TF_exp_only \
                   --add_predict_TF $add_predict_tf --add_YRI $add_YRI --population $population \
                   --TF_exp_type $TF_exp_type --add_gm12878 $add_gm12878 --new_batch $new_batch \
                   --batch_mode $new_batch_random --other_info ${17}
else
    echo 'Classification'
    Rscript3 $project_dir/s_classification_for_one_gene.R --batch_name $1 --test FALSE --add_histone TRUE --add_miRNA TRUE --model $2 --gene $3
fi














