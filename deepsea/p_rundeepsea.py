#python rundeepsea.py infilename outdir
#The output files will be in outdir.

from subprocess import *
from tempfile import mkdtemp
import sys

import os
home_dir = os.path.expanduser('~')
lib_dir = '%s/python/' % home_dir
import sys
sys.path.insert(0, lib_dir)
sys.path.insert(0, '%s/expression_var/python/' % home_dir)
import pandas as pd
import p_mymodule as my
from p_project_metadata import *

import argparse
parser = argparse.ArgumentParser(description='Extract the deepsea predictions of one tf.Deepsea prediction is sample based. The output of this script is TF based.')
if __doc__ is None:
    parser.add_argument('--vcf_file', help = "vcf file of the chromose", default =None)
    parser.add_argument('--out_dir', help = "Output dir", default =None)
    args = parser.parse_args()
    vcf_file = args.vcf_file
    outdir = args.out_dir
else:
    vcf_file = './examples/deepsea/example.vcf'
    outdir = 'outdir'
    
cpoutdir=True
check_call(['mkdir','-p',outdir])

peak_file_df_rmdup = f_get_peak_file_df_rmdup(project_dir, 'processed')



if 'vcf' in vcf_file:
    
    for loc_tf in peak_file_df_rmdup.tf:
        final_file = '%s/%s.out.evalue' % (outdir, loc_tf)
        if os.path.isfile(final_file):
            logging.info('Skip %s' % loc_tf)
            continue

        try:
            tempdir = mkdtemp()
            tmp_dir = tempdir
            peak_file = peak_file_df_rmdup.ix[loc_tf, 'file_path']
            #vcf_file = '%s/deepsea/tests/data/chr22.merge.head.vcf.gz'%(project_dir)
            deepsea_tf = peak_file_df_rmdup.ix[loc_tf, 'deepsea_tf']
            
            print "Successfully copied input to working directory " + tempdir 
            try:
                my.f_shell_cmd("python2.7 p_generate_peak_fastq.py --vcf_file %s --peak_file %s --tmp_dir %s" % (vcf_file, peak_file, tmp_dir))
            except:
                raise Exception('Vcf format error.')
            #retrieve 1100bp instead of 1000bp for supporting deletion variants (<100bp) 
            check_call(["python2.7 p_fasta2input.py --fasta_file %s/infile.vcf.wt1100.fasta"% tmp_dir],shell=True) 
            print "Successfully converted to input format"
        
        
            check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.vcf.wt1100.fasta.ref.h5"],shell=True)
            check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.vcf.mut1100.fasta.ref.h5"],shell=True)
            print "Finished running DeepSEA. Now prepare output files..."
            my.f_shell_cmd("python p_h5ToOutput.py %s/infile.vcf %s/infile.vcf.wt1100.fasta.ref.h5.pred.h5  %s/infile.vcf.mut1100.fasta.ref.h5.pred.h5  %s" %(tmp_dir, tmp_dir, tmp_dir, deepsea_tf)) 
            
            if cpoutdir:
                my.f_ensure_make_dir(outdir)
                my.f_shell_cmd("cp %s/infile.vcf.out.ref %s/%s.out.ref" % (tmp_dir, outdir, loc_tf))
                my.f_shell_cmd("cp %s/infile.vcf.out.diff %s/%s.out.diff" % (tmp_dir, outdir, loc_tf))
                my.f_shell_cmd("cp %s/infile.vcf.out.evalue %s/%s.out.evalue" % (tmp_dir, outdir, loc_tf))
        
        except:
            print 'Skip %s' % loc_tf
        print "Finished creating output file. Now clean up..."
        call(['rm',tempdir,'-r'])
        print "Everything done."





