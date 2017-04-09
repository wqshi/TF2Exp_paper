from optparse import OptionParser
import os
import math
import sys
sys.path.append("/homed/home/shi/python/")
from time import gmtime, strftime
import p_mymodule as my

project_dir = os.path.expanduser( '~/projects/expression_var/sailfish/data/')
script_dir = os.path.expanduser('~/projects/expression_var/sailfish/')


#my.f_shell_cmd('scp /homed/home/shi/sailfish/quantify_RNA_node_dir.py $CLUST:/home/shi/expression_var/sailfish/')

def main():
    parser=OptionParser()

    if __doc__ is None:
        parser.add_option('--out',dest='out',help='Out',default='%s/qsub/'%project_dir)
        opts,args=parser.parse_args()
        out_dir = opts.out
        node_dir="/state/partition1/shi/tmp_depth/%s/" % my.f_shell_cmd('echo $JOB_ID', quiet = True).replace('\n', '')
    else:
        out_dir = '%s/qsub/'%project_dir
        node_dir = out_dir + '/node/'

    my.f_ensure_make_dir(node_dir)
    FQ_dir='%s/fastq/' % project_dir
    geuvadis_meta='%s/metaData/E-GEUV-1.sdrf.txt' % project_dir
    our_study='%s/metaData/our_sample.list1' % project_dir
    metadata='%s/metadata' % project_dir
    
    
    our_people=set()
    gender={}
    for line in open(our_study,'r').readlines():
        our_people.add(line.strip().split('\t')[0])
        items=line.strip().split('\t')
        person=items[0]
        person_gender=items[3]
        if person not in gender.keys():
            gender[person]=person_gender
      
    geu1=set()
    for line in open(geuvadis_meta,'r').readlines():
        items=line.strip().split('\t')
        geu1.add(items[0])

    of_interest=geu1.intersection(our_people)
    print of_interest
    print len(of_interest)

    person_to_fq={}
    for line in open(geuvadis_meta,'r').readlines():
        items=line.strip().split('\t')
        person=items[0]
        if person not in of_interest:
            continue
        if person not in person_to_fq.keys():
            person_to_fq[person]=set()
        curr_fq=items[28]
        person_to_fq[person].add(FQ_dir+os.path.basename(curr_fq))
        print items
    
    print person_to_fq
    metadata_file=open(metadata,'w')
    for person in person_to_fq.keys(): 
        
        out_curr= node_dir + person+'.sailfish'
        metadata_file.write(person+'\t'+','.join(person_to_fq[person])+'\t'+out_curr+'\n')
        #And run sailfish
        cur_gender=gender[person]
        sailfish_idx='%s/Transcriptome/gencode.v19.annotation.PC.lincRNA.gtf.splicedExon.N'% project_dir +cur_gender+'.fa.dedup.fa_IDX_sailfish'
        #cmd_module='module load sailfish/0.6.3'
        library_type='"T=PE:O=><"' #T=PE:O=><:S=SA
        fastqs=list(person_to_fq[person])
        cmds=[]
        cmds.append('#!/usr/bin/env bash')
        cmds.append('cp %s %s' % ( ' '.join(fastqs), node_dir )  )
        loc_fastqs = [ os.path.join(node_dir, os.path.basename(fastq_file)) for fastq_file in fastqs ]

        #cmds.append(cmd_module)
        sailfish_cmd='sailfish quant -i '+sailfish_idx+' -l '+library_type+' -1 <(gunzip -c '+loc_fastqs[0]+') -2 <(gunzip -c '+loc_fastqs[1]+') -o '+out_curr+' -f'
        cmds.append(sailfish_cmd)
        cmds.append('cd '+out_curr)
        #cmds.append('module load java/latest')
        gtf='%s/GENCODE_v19_2014-06-03/gencode.v19.annotation.PC.lincRNA.gtf' % project_dir
        cmds.append('%s/TranscriptsToGenes.sh --gtf-file '% script_dir +gtf+' --exp-file '+out_curr+'/quant.sf'+' --res-file '+person+'quant.gene_level.sf')
        cmds.append('mv '+out_curr+'/quant.sf'+' '+out_curr+'/'+person+'quant.sf')
        #cmds.append('/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/src/TranscriptsToGenes.sh --gtf-file '+gtf+' --exp-file '+out_curr+'/quant_bias_corrected.sf'+' --res-file '+person+'quant_bias_corrected.gene_level.sf')
        cmds.append('mv %s %s/' %(out_curr, out_dir ))
        cmds.append('\rm %s' %(out_curr))
        print '\n'.join(cmds)
        qsub_a_command('qqqq'.join(cmds),out_curr+'_script.sh','qqqq','10G')

def qsub_a_command(cmd,shell_script_name,split_string=',',memory_number='20G'):
    f=open(shell_script_name,'w')
    print shell_script_name
    cmds=cmd.split(split_string)
    for i in range(len(cmds)):
        f.write(cmds[i]+'\n')
    f.close()
    os.system('chmod 711 '+shell_script_name)
    #os.system("qsub -V -q shi.q -l mem_free="+memory_number+" -l h_vmem="+memory_number+" -l h_rt=20:00:00 -o "+shell_script_name+'.o'+' -e '+shell_script_name+'.e'+' '+shell_script_name)



main()
