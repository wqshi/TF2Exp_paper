pop_name=$1
head_dir=~/expression_var/data/raw_data/pop/$pop_name'_dir'/
#RNAdata=~/expression_var/data/raw_data/pop/$pop_name'_dir'/node
RNAdata=/state/partition1/shi/pop/$pop_name'_dir'/
RNAdataSNP=${RNAdata}SNPdata
RNAdataTrans=${RNAdata}Transcriptome
mkdir -p ${RNAdata}
mkdir ${RNAdataTrans}
mkdir ${RNAdataSNP}

#========================================
#1. Make a genome with Ns instead of SNPs
#======================================== 
#- vcf with N as alternative
#high coverage seq

vcf=~/expression_var/data/raw_data/pop/$pop_name'_dir'/chr$pop_name.vcf.gz
zcat ${vcf} | grep -v '##' | awk '{print $1"\t"$2"\t"$3"\t"$4"\tN\t"$6"\t"$7"\t"$8"\t"$9"\t1|1"}' > ${RNAdataSNP}/SNPs_for_N.vcf
#Remove duplicate SNPs from the large vcf file
cat ${RNAdataSNP}/SNPs_for_N.vcf | awk '!_[$3]++' > ${RNAdataSNP}/SNPs_for_N.dedup.vcf
cat ${RNAdataSNP}/SNPs_for_N.dedup.vcf | sed 's/FORMAT\t1|1/FORMAT\tNperson/g' > ${RNAdataSNP}/SNPs_for_N.dedup.indName.vcf
rm ${RNAdataSNP}/SNPs_for_N.dedup.vcf
#Add Ns to genome - add to the male genome
p=python2.7
code=~/expression_var/python/addSNPsToFa.py
vcf=${RNAdataSNP}/SNPs_for_N.dedup.indName.vcf
fadir=~/projects/wgs/hg19ByChrom/maleByChrom/
dict=~/expression_var/data/raw_data/pop/hg19.fa.fai
outpref=${RNAdataTrans}/N_hg19_ENCODE
indiv=Nperson
${p} ${code} -v ${vcf} --unphased '' ${fadir} ${dict} ${outpref} ${indiv}
#DONE! We have a genome with Ns.

#We'll have 1 universal male genome, and we'll get female/male specific transcriptomes by subsetting gtf file below.
Ngenome=${RNAdataTrans}/Ngenome.hg19_ENCODE.male.fa
cp ${outpref}.paternal.fa ${Ngenome}
rm ${outpref}.paternal.fa
rm ${outpref}.maternal.fa
#==========================================

#======================================================================================
#2. gffread to make fasta genome (maternal=female and paternal=male) into transcriptome
#======================================================================================
cp ~/expression_var/data/raw_data/pop/gencode.v19.annotation.gtf.gz $RNAdata
gencode_gtf=${RNAdata}/gencode.v19.annotation.PC.lincRNA.gtf
gencode_gtf_gz=${RNAdata}/gencode.v19.annotation.gtf.gz
zcat ${gencode_gtf_gz} | grep "protein_coding\|lincRNA" > ${gencode_gtf}
#Female transcriptome
gtf_female=${RNAdata}/gencode.v19.annotation.female.PC.lincRNA.gtf
cat ${gencode_gtf} | grep -v chrY > ${gtf_female}
out_female=${RNAdataTrans}/gencode.v19.annotation.PC.lincRNA.gtf.splicedCDS.Nfemale.fa.out
splicedCDS_female=${RNAdataTrans}/gencode.v19.annotation.PC.lincRNA.gtf.splicedCDS.Nfemale.fa
splicedExon_female=${RNAdataTrans}/gencode.v19.annotation.PC.lincRNA.gtf.splicedExon.Nfemale.fa

cuff_dir=~/packages/cufflinks-2.2.0.Linux_x86_64/
$cuff_dir/gffread ${gtf_female} -g ${Ngenome} -s ${dict} -x ${splicedCDS_female} -w ${splicedExon_female} -o ${out_female}
${p} ~/expression_var/sailfish/FilterByName.py --input ${splicedExon_female} --output ${splicedExon_female}.dedup.fa

#Male transcriptome
out_male=${RNAdataTrans}/gencode.v19.annotation.PC.lincRNA.gtf.splicedCDS.Nmale.fa.out
splicedCDS_male=${RNAdataTrans}/gencode.v19.annotation.PC.lincRNA.gtf.splicedCDS.Nmale.fa
splicedExon_male=${RNAdataTrans}/gencode.v19.annotation.PC.lincRNA.gtf.splicedExon.Nmale.fa

$cuff_dir/gffread ${gencode_gtf} -g ${Ngenome} -s ${dict} -x ${splicedCDS_male} -w ${splicedExon_male} -o ${out_male}
${p} ~/expression_var/sailfish/FilterByName.py --input ${splicedExon_male} --output ${splicedExon_male}.dedup.fa

#======================
#3. Index transcriptome
#======================
splicedExon_female=${RNAdataTrans}/gencode.v19.annotation.PC.lincRNA.gtf.splicedExon.Nfemale.fa
splicedExon_male=${RNAdataTrans}/gencode.v19.annotation.PC.lincRNA.gtf.splicedExon.Nmale.fa

#Exons
#=====
#Female                           
splicedExon_female_dedup=${splicedExon_female}.dedup.fa                                                                                                                      
out_idx=${splicedExon_female_dedup}_IDX_sailfish
script_idx=${out_idx}_script.sh
#echo "module load sailfish/0.6.3"> ${script_idx}
sailfish_dir='~/packages/Sailfish-0.6.3-Linux_x86-64/bin/'
echo "$sailfish_dir/sailfish index -t ${splicedExon_female_dedup} -o ${out_idx} -k 24" >> ${script_idx}
chmod 711 ${script_idx}
sh ${script_idx}
#qsub -l h_vmem=20G -o ${script_idx}.o -e ${script_idx}.e ${script_idx}
cp -r $out_idx $head_dir

#Male                                                                                                                                                                       
splicedExon_male_dedup=${splicedExon_male}.dedup.fa
out_idx=${splicedExon_male_dedup}_IDX_sailfish
script_idx=${out_idx}_script.sh
#echo "module load sailfish/0.6.3"> ${script_idx}
echo "$sailfish_dir/sailfish index -t ${splicedExon_male_dedup} -o ${out_idx} -k 24 -f" >> ${script_idx}
chmod 711 ${script_idx}
sh ${script_idx}
#qsub -l h_vmem=20G -o ${script_idx}.o -e ${script_idx}.e ${script_idx}

cp -r $out_idx $head_dir
rm -r $RNAdata
