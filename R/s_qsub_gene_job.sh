script_dir='/home/shi/expression_var/R/'
#batch_name='462samples_sailfish_quantile'
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
#sh $HOME/expression_var/R//s_rsync_regression_data.sh $batch_name
echo $PWD
#ssh shi@clustdell.cmmt.ubc.ca /home/shi/s_qsub_shell.sh $model-$batch_name-shi-24G $script_dir/s_run_single.sh $batch_name $model

gene_dir=$HOME/expression_var/R/data/$batch_name/rnaseq/$chr_str/
cd $gene_dir
#for gene_file in $(ls *.txt)
for gene_file in $(find *.txt -size +0)
do
    gene_name=${gene_file%.txt}

    #if grep -q ${gene_name%.*} $HOME/expression_var/R//data/r_results/max_perf.txt; then
        sh /home/shi/s_qsub_shell.sh $model-$chr_str-$batch_name-${16}-$gene_name-shi-small $script_dir/s_run_single_gene.sh $batch_name $test_flag $model $chr_str $gene_name $add_miRNA $add_TF_exp $add_permutation $add_TF_exp_only $add_predict_tf $add_YRI ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19}
  
        if [ $test_flag == 'TRUE' ];then
            break
        fi

    #else
    #    echo "skip $gene_name"
    #fi
done

sh /home/shi/s_qsub_python.sh summary-$chr_str-$batch_name-${16}-$gene_name-shi-small $script_dir/p_detect_job_finished_run_new.py --batch_name $batch_name --last_gene $gene_name --target_mode ${16} --chr_str $chr_str --add_penalty $add_permutation --other_info ${17} --new_batch ${15}






