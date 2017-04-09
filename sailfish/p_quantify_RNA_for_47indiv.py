
import os
import math
import sys
sys.path.append("/homed/home/shi/python/")
from time import gmtime, strftime
import p_mymodule as my

project_dir = os.path.expanduser( '~/projects/expression_var/sailfish/data/')
script_dir = os.path.expanduser('~/projects/expression_var/sailfish/')
import argparse

#my.f_shell_cmd('scp /homed/home/shi/sailfish/quantify_RNA_node_dir.py $CLUST:/home/shi/expression_var/sailfish/')
parser = argparse.ArgumentParser(description='Extract the deepsea predictions of one tf.Deepsea prediction is sample based. The output of this script is TF based.')

def main():

    if __doc__ is None:
        parser.add_argument('--out_dir',help='Out',default='%s/qsub_47samples/'%project_dir)
        parser.add_argument('--test_flag',help='Test flag',default='T')
        opts = parser.parse_args()
        out_dir = opts.out_dir
        test_flag = (opts.test_flag == 'T')
        node_dir="/state/partition1/shi/tmp_depth/%s/" % my.f_shell_cmd('echo $JOB_ID', quiet = True).replace('\n', '')
    else:
        out_dir = '%s/qsub_47samples/'%project_dir
        node_dir = out_dir + '/node/'
        test_flag = True

    my.f_ensure_make_dir(out_dir)
    FQ_dir='%s/fastq47indiv/' % project_dir
    geuvadis_meta='%s/metaData/E-MTAB-3656.sdrf.txt' % project_dir
    our_study='%s/metaData/our_sample.list' % project_dir
    metadata='%s/metadata' % project_dir
    
    #import ipdb; ipdb.set_trace()
    our_people=set()
    gender={}
    pop={}
    for line in open(our_study,'r').readlines():
        our_people.add(line.strip().split('\t')[0])
        items=line.strip().split('\t')
        person=items[0]
        person_gender=items[3]
        if person not in gender.keys():
            gender[person]=person_gender
        if person not in pop.keys():
            pop[person]=items[1]
            
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

        person_to_fq[person].add(FQ_dir+'/%s_end1.fastq.gz'%person )
        person_to_fq[person].add(FQ_dir+'/%s_end2.fastq.gz'%person )
        #print items

    #import ipdb; ipdb.set_trace()
    print person_to_fq
    metadata_file=open(metadata,'w')
    for person in person_to_fq.keys(): 
        
        out_curr= node_dir + person+'.sailfish/'
        metadata_file.write(person+'\t'+','.join(person_to_fq[person])+'\t'+out_curr+'\n')
        #And run sailfish
        cur_gender=gender[person]
        cur_pop = pop[person]
        #sailfish_idx='%s/Transcriptome/gencode.v19.annotation.PC.lincRNA.gtf.splicedExon.N'% project_dir +cur_gender+'.fa.dedup.fa_IDX_sailfish'
        index_dir = '~/expression_var/data/raw_data/pop/%s_dir' % cur_pop
        sailfish_idx='%s/gencode.v19.annotation.PC.lincRNA.gtf.splicedExon.N'% index_dir +cur_gender+'.fa.dedup.fa_IDX_sailfish'
        #cmd_module='module load sailfish/0.6.3'
        library_type='"T=PE:O=><"' #T=PE:O=><:S=SA
        fastqs=list(person_to_fq[person])

        #If the output is there, don't lanch the jobs again.
        final_out_file = '%s/%s.sailfish/%squant.gene_level.sf' % ( out_dir, person, person)
        if os.path.isfile(final_out_file):
            print 'Got the results of %s' % person
            continue
        else:
            print 'Sailfish %s' % person
            my.f_remove_dir( '%s/%s.sailfish' % ( out_dir, person))

        if not os.path.isfile(fastqs[0]):
            print 'Missing person %s: %s'% (person, fastqs[0])
            if not os.path.isfile(fastqs[1]):
                print 'Missing person %s: %s'% (person, fastqs[1])
                continue
            continue
        cmds=[]
        cmds.append('#!/usr/bin/env bash')
        cmds.append('mkdir -p %s' % out_curr)
        cmds.append('cp -u %s %s' % ( ' '.join(fastqs), out_curr )  )
        loc_fastqs = [ os.path.join(out_curr, os.path.basename(fastq_file)) for fastq_file in fastqs ]
        loc_fastqs.sort()
        #import ipdb; ipdb.set_trace()
        #cmds.append(cmd_module)
        sailfish_exe='~/packages/Sailfish-0.6.3-Linux_x86-64/bin/sailfish'
        sailfish_cmd=sailfish_exe +' quant -i '+sailfish_idx+' -l '+library_type+' -1 <(gunzip -c '+loc_fastqs[0]+') -2 <(gunzip -c '+loc_fastqs[1]+') -o '+out_curr+' -f'
        cmds.append(sailfish_cmd)
        cmds.append('cd '+out_curr)
        #cmds.append('module load java/latest')
        gtf='%s/GENCODE_v19_2014-06-03/gencode.v19.annotation.PC.lincRNA.gtf' % project_dir
        cmds.append('%s/TranscriptsToGenes.sh --gtf-file '% script_dir +gtf+' --exp-file '+out_curr+'/quant.sf'+' --res-file '+person+'quant.gene_level.sf')
        cmds.append('mv '+out_curr+'/quant.sf'+' '+out_curr+'/'+person+'quant.sf')
        cmds.append('rm %s/*.fastq.gz' %(out_curr))
        cmds.append('rm %s/reads.*' %(out_curr))
        cmds.append('mv %s %s/' %(out_curr, out_dir ))
        cmds.append('rm -r %s' %(out_curr))
        print '\n'.join(cmds)
        if test_flag == False:
            qsub_a_command('qqqq'.join(cmds),out_dir + person +'_script.sh','qqqq','10G')
        else:
            print 'Skip %s' % person

def qsub_a_command(cmd,shell_script_name,split_string=',',memory_number='20G'):
    f=open(shell_script_name,'w')
    print shell_script_name
    cmds=cmd.split(split_string)
    for i in range(len(cmds)):
        f.write(cmds[i]+'\n')
    f.close()
    os.system('chmod 711 '+shell_script_name)
    #import ipdb; ipdb.set_trace()
    if my.f_get_server_name()=='loire':
        os.system('sh %s'%shell_script_name)
    else:
        os.system("qsub -V -q shi.q -l mem_free="+memory_number+" -l h_vmem="+memory_number+" -l h_rt=20:00:00 -o "+shell_script_name+'.o'+' -e '+shell_script_name+'.e'+' '+shell_script_name)

main()

#p quantify_RNA_node_dir.py qsub_445samples #Run in the clustdell head dir.













