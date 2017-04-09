
pop_list=(CEU  FIN  GBR  TSI  YRI)
pop_list=(FIN)
script_dir='/home/shi/expression_var/sailfish/'
sh ../python/s_rsnyc_to_cluster.sh 

for loc_pop in ${pop_list[@]}
do
  echo $loc_pop  
  ssh shi@clustdell.cmmt.ubc.ca /home/shi/s_qsub_shell.sh index-shi-24G-$loc_pop  $script_dir/RNA_mapping.sh $loc_pop
  #break
done
