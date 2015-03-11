'''
Created on July 18, 2011
@author: dazhe.meng

Updated imputation script to run mach
Newest update: add new mode: sp_prepare_ref, sp_prepare_samp, for sparse sequencing 
'''

import subprocess
import os, sys
import cPickle

from ConfigParser import ConfigParser
cfgs = ConfigParser()
cfg_details="""
Dazhe Configuration File Details:
---------------------------------------------
The .dazhe_config file should be placed in the home directory
It is similar in format to ini files, here is a sample (with required params!)
[PATH]
home_dir = /home/cmb-03/mn/dazhemen/
snpdata_dir = %(home_dir)s/SNPdata/
jobscript_dir = %(home_dir)s/jobscripts/
[IMPUTE]
reference_data = 
sample_data = 250K_m75_chr%s.csv
min_maf = 0.05
max_mis = 0.30
[MACH]
script_dir = %(home_dir)s/Impute/
program_dir = %(home_dir)s/Programs/MACH/
output_dir = %(home_dir)s/MACH/
chunk_size = 5000
num_iter = 30
# The processed dir stores filtered SNP data as well as parameter information
processed_dir = %(home_dir)s/MACH/
# The data identifier uniquely specifies which dataset is used in preprocessing
data_identifier = swe180_250K_%(chunk_size)s
[USCCLUSTER]
mem_required = 3000mb
wall_time = 24:00:00
"""

#cfgs.read(os.path.expanduser('~/.dazhe_config'))
#sys.path.append(cfgs.get('PATH','pylibrary_dir'))

from argparse import ArgumentParser
parser = ArgumentParser(description = "Run imputation on cluster", epilog="", prefix_chars="-")
parser.add_argument('mode', choices=['preprocess','impute','mle','read','error','sp_prepare_ref','sp_prepare_samp','readdose'], default='impute')
parser.add_argument("-c", "--config_file", help="Name of the configuration file to use", default="~/.dazhe_config")
parser.add_argument("-r", "--reference_data", help="Name of the data with more snps", default="")
parser.add_argument("-i", "--sample_data", help="Name of the data with more accessions", default="")
parser.add_argument("-n", "--impute_job_name", required = True, help="name of this imputation job, doubles as output dir and temporary file name")
parser.add_argument("-a", "--sparse_acc", help="used in sp_prepare_ref to set which acc to hide", default="")
parser.add_argument("-d", "--data_identifier", help="Data identifier", default="")
parser.add_argument("--merged_chr", action='store_true', help="specifies that data is not separated into chromosomes")
parser.add_argument("--test_run", type=int, nargs='?', help="specifies the number of jobs to create for a testing run", default=0, const=1)
parser.add_argument("--mem_req", help="(USC Cluster) Modify the memory requirement", default="")
args = parser.parse_args()

cfgs.read(os.path.expanduser(args.config_file))
cfgs.set('MACH','home_dir', cfgs.get('PATH','home_dir')) # copy the path information!
cfgs.set('MACH','script_home_dir', cfgs.get('PATH','script_home_dir')) # copy extra path information for GMI!
if not args.data_identifier:
    cfgs.set('MACH','data_identifier', args.impute_job_name)
if not os.path.isdir(cfgs.get('MACH','output_dir')):
    os.mkdir(cfgs.get('MACH','output_dir'))
if not os.path.isdir(cfgs.get('MACH','output_dir')+cfgs.get('MACH','data_identifier')):
    os.mkdir(cfgs.get('MACH','output_dir')+cfgs.get('MACH','data_identifier'))
if not os.path.isdir(cfgs.get('MACH','output_dir')+cfgs.get('MACH','data_identifier')+'/sep/'):
    os.mkdir(cfgs.get('MACH','output_dir')+cfgs.get('MACH','data_identifier')+'/sep/')
if args.mem_req != "":
    cfgs.set('USCCLUSTER','mem_required', args.mem_req)
if args.sample_data != "":
    cfgs.set('IMPUTE', 'sample_data', args.sample_data)
if args.reference_data != "":
    cfgs.set('IMPUTE', 'reference_data', args.reference_data)
# writes a new config file
new_cfg_filename = '%s%s.cfg'%(cfgs.get('MACH','output_dir'), args.impute_job_name)
if os.path.isfile(new_cfg_filename):
    print "Config file existing, not overwriting old config"
else:
    with open(new_cfg_filename, 'w') as newconfigfile:
        cfgs.write(newconfigfile)

import socket
hostname = socket.gethostname()
jf_header = ""
if hostname == 'hpc-cmb.usc.edu':  #USC
    host = "usc"
    jf_header += "#!/bin/bash\n"
    jf_header += "#PBS -l walltime=%s\n"%cfgs.get('USCCLUSTER','wall_time')
    jf_header += "#PBS -l mem=%s\n"%cfgs.get('USCCLUSTER','mem_required')
    jf_header += "#PBS -q cmb\n"
    jf_header += "#PBS -k eo\n"
    jf_header += "export PYTHONPATH='/home/cmb-01/mn/dazhemen/PyLibrary/:/home/cmb-01/mn/dazhemen/Biopython/'\n"
    jf_header += "export PATH='/usr/usc/python/2.6.1/bin/':$PATH\n"
elif hostname == 'mgmt01': #GMI
    host = "gmi"
    jf_header += "#!/bin/sh\n"
    jf_header += "source /etc/modules-env.sh\n"
    jf_header += "module load matplotlib/1.0.0\n"
    jf_header += "module load numpy/MKL/1.5.1\n"
    jf_header += "export PYTHONPATH=$PYTHONPATH:'/home/GMI/dazhe.meng/PyLibrary/'\n"

# convenience constants so as not to refer to the configuration object all the time
script_dir = cfgs.get('MACH','script_dir')
output_dir = cfgs.get('MACH','output_dir')+cfgs.get('MACH', 'data_identifier')+"/"

if args.mode == 'preprocess':
    # The preprocessing step should remove problematic snps as well as separate the chromosome into chunks
    if 0 and os.path.isfile(output_dir+"chunk_info.pickled"):
        print "Preprocessing for this dataset already finished!"
    else:
        for chri in xrange(1,6):
            jn = 'mach_preprocess_%s_%s'%(cfgs.get('MACH', 'data_identifier'), chri)
            jfn = cfgs.get('PATH','jobscript_dir')+jn+".sh"
            jf = open(jfn, "w")
            jf.write(jf_header)
            jf.write('#$ -pe threads 2\n')
            # Fix this line later!
            if args.sparse_acc:
                jf.write("python %sRun_MACH.py preprocess -n %s -c %s -a %s\n"%(script_dir, new_cfg_filename, chri, args.sparse_acc))
            else:
                jf.write("python %sRun_MACH.py preprocess -n %s -c %s\n"%(script_dir, new_cfg_filename, chri))
            jf.close()
            subprocess.call("qsub %s"%jfn, shell=True)
elif args.mode == 'impute':
    # The impute step doubles as parameter estimator
    if 0 and not os.path.isfile(output_dir+"chunk_info.pickled"):
        # the chunk info file is a pickled list of chunk strings <chr>_<chunks_ind>
        print "Preprocessing not found, please use preprocess first"
    else:
        chunks = cPickle.load(open(output_dir+"chunk_info.pickled",'rb'))
        for chunk in chunks:
            if os.path.isfile(output_dir+'sep/%s.dose.gz'%chunk):
                print "done: %s"%chunk
                continue
            jn = 'mach_impute_%s_%s'%(cfgs.get('MACH', 'data_identifier'), chunk)
            jfn = cfgs.get('PATH','jobscript_dir')+jn+".sh"
            jf = open(jfn, "w")
            jf.write(jf_header)
            jf.write('#$ -pe threads 2\n')
            jf.write('#$ -q q.norm@blade*\n')
            jf.write("python %sRun_MACH.py impute -n %s -f %s\n"%(script_dir, new_cfg_filename, chunk))
            jf.close()
            subprocess.call("qsub %s"%jfn, shell=True)
elif args.mode == 'mle':
    # Fix later
    pass
elif args.mode == 'read':
    for chri in xrange(1,6):
        jn = 'mach_read_%s_%s'%(cfgs.get('MACH', 'data_identifier'), chri)
        jfn = cfgs.get('PATH','jobscript_dir')+jn+".sh"
        jf = open(jfn, "w")
        jf.write(jf_header)
        # Fix this line later!
        jf.write("python %sRun_MACH.py read -n %s -c %s\n"%(script_dir, new_cfg_filename, chri))
        jf.close()
        subprocess.call("qsub %s"%jfn, shell=True)
elif args.mode == 'readdose':
    for chri in xrange(1,6):
        jn = 'mach_readdose_%s_%s'%(cfgs.get('MACH', 'data_identifier'), chri)
        jfn = cfgs.get('PATH','jobscript_dir')+jn+".sh"
        jf = open(jfn, "w")
        jf.write(jf_header)
        jf.write('#$ -pe threads 2\n')
        jf.write('#$ -q q.norm@blade*\n')
        # Fix this line later!
        jf.write("python %sRun_MACH.py readdose -n %s -c %s\n"%(script_dir, new_cfg_filename, chri))
        jf.close()
        subprocess.call("qsub %s"%jfn, shell=True)

