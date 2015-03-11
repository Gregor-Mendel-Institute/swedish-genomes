"""
Script that defines job on a single processor
The environment for this script should have been setup by the jobscript
"""

import os, sys
import subprocess
import cPickle
from ImputeData import ImputeData
from ConfigParser import ConfigParser
cfgs = ConfigParser()

NUM_CHR = 5

from argparse import ArgumentParser
parser = ArgumentParser(description = "Run imputation on cluster", epilog="", prefix_chars="-")
parser.add_argument('mode', choices=['preprocess','impute','mle','read','error','sp_prepare_ref','sp_prepare_samp','readdose'], default='impute')
parser.add_argument("-n", "--config_file", help="Name of the configuration file to use", default="~/.dazhe_config")
parser.add_argument("-f", "--file_prefix", help="(For impute, mle) Name of the file(s) to operate on")
parser.add_argument("-a", "--sparse_acc", help="Used in sparse modes", default="")
parser.add_argument("-c", "--chr", type=int, help="(For preprocess, read) Which chromosome to operate on")
args = parser.parse_args()

result_list = []
cfgs.read(os.path.expanduser(args.config_file))
script_dir = cfgs.get('MACH','script_dir')
snpdata_dir = cfgs.get('PATH','snpdata_dir')
output_dir = cfgs.get('MACH','output_dir')+cfgs.get('MACH', 'data_identifier')+'/'
snp_bin_size = cfgs.getint('MACH','chunk_size')
sep_output_dir = output_dir + '/sep/'

if args.mode == 'preprocess':
    print "I'm running on chr %s, commencing preprocessing"%(args.chr)
    I = ImputeData(snpdata_dir+cfgs.get('IMPUTE','reference_data')%args.chr) # The reference data should be set no matter what!
    if args.sparse_acc != "":
        print "Removing acc"
        I.genotype.remove_acc(args.sparse_acc)
    if cfgs.get('IMPUTE','sample_data')=='self': # Might have to fix later
        I.filter_snps(max_mis = cfgs.getfloat('IMPUTE','max_mis'), min_maf=cfgs.getfloat('IMPUTE','min_maf'), max_allele=2)
        I.Chromosome_divideup(snp_bin_size, snp_bin_size/10, snp_bin_size/5*2)
        I.MACHoutput(sep_output_dir+str(args.chr), 'Single')
        I.genotype.output(output_dir+"filtered_ref_chr%s.csv"%args.chr) # also output the filtered genotype file for future use
        I.Output_chunks(output_dir+"chunk_chr%s"%args.chr)       
    else:
        I.filter_snps(max_mis = cfgs.getfloat('IMPUTE','max_mis'), min_maf=cfgs.getfloat('IMPUTE','min_maf'))
        I.add_imputing_accs(snpdata_dir+cfgs.get('IMPUTE','sample_data')%args.chr)
        I.filter_snps(max_allele=2)
        I.Chromosome_divideup(snp_bin_size, snp_bin_size/10, snp_bin_size/5*2)
        I.separate_data()
        I.genotype.output(output_dir+"filtered_ref_chr%s.csv"%args.chr) # also output the filtered genotype file for future use
        I.target_genotype.output(output_dir+"filtered_samp_chr%s.csv"%args.chr) # also output the filtered genotype file for future use
        I.MACHoutput(sep_output_dir+str(args.chr), 'Double')
        I.Output_chunks(output_dir+"chunk_chr%s"%args.chr)
    # Check if all chunk file are there, and then combine them if they do, an ugly way to avoid MPI, and should work if things finish at different times
    num_chunks = [0]*NUM_CHR
    if 0 in num_chunks:
        print 'yeah'
    for chri in xrange(1,NUM_CHR+1):
        if os.path.isfile(output_dir+"chunk_chr%s"%chri):
            if os.stat(output_dir+"chunk_chr%s"%args.chr).st_size>0: # is not empty
                num_chunks[chri-1] = int(open(output_dir+"chunk_chr%s"%chri).readline().strip())
    if 0 not in num_chunks: # we find all files!
        print num_chunks
        chunks = []
        for chri, nchunk in enumerate(num_chunks):
            for ichunk in xrange(nchunk):
                chunks.append('%s_%s'%(chri+1,ichunk))
            # cleanup
            subprocess.call('rm '+output_dir+"chunk_chr%s"%(chri+1), shell=True)
        cPickle.dump(chunks, open(output_dir+"chunk_info.pickled","wb"))
elif args.mode == 'sp_prepare_ref':
    print "I'm running on chr %s, commencing sp ref preprocessing"%(args.chr)
    I = ImputeData(snpdata_dir+cfgs.get('IMPUTE','reference_data')%args.chr) # The reference data should be set no matter what!
    I.genotype_data.remove_acc(args.sparse_acc)
    I.Chromosome_divideup(snp_bin_size, snp_bin_size/10, snp_bin_size/5*2)
    I.MACHoutput(sep_output_dir+str(args.chr), 'Single')
    I.Output_chunks(output_dir+"chunk_chr%s"%args.chr)       
    num_chunks = [0]*NUM_CHR
    if 0 in num_chunks:
        print 'yeah'
    for chri in xrange(1,NUM_CHR+1):
        if os.path.isfile(output_dir+"chunk_chr%s"%chri):
            if os.stat(output_dir+"chunk_chr%s"%args.chr).st_size>0: # is not empty
                num_chunks[chri-1] = int(open(output_dir+"chunk_chr%s"%chri).readline().strip())
    if 0 not in num_chunks: # we find all files!
        print num_chunks
        chunks = []
        for chri, nchunk in enumerate(num_chunks):
            for ichunk in xrange(nchunk):
                chunks.append('%s_%s'%(chri+1,ichunk))
            # cleanup
            subprocess.call('rm '+output_dir+"chunk_chr%s"%(chri+1), shell=True)
        cPickle.dump(chunks, open(output_dir+"chunk_info.pickled","wb")) 
elif args.mode == 'impute':
    cmd = ["%smach1"%cfgs.get('MACH','program_dir'), "--dosage","--greedy", "-r %s"%cfgs.getint('MACH','num_iter')]
    cmd += ["-d", sep_output_dir+args.file_prefix+".dat"]
    cmd += ["-p", sep_output_dir+args.file_prefix+".ped"]
    if cfgs.get('IMPUTE','sample_data')!='self':
        cmd += ["-s", sep_output_dir+args.file_prefix+".snps"]
        cmd += ["-h", sep_output_dir+args.file_prefix+".haplos"]
    cmd += ["--prefix", sep_output_dir+args.file_prefix]
    retcode = subprocess.call(cmd)             
elif args.mode == 'mle': # note that this should never be run when there is no separate target genotypes!
    cmd = ["%smach1"%cfgs.get('MACH','program_dir'), "--greedy", "--mle"]
    cmd += ["-d", sep_output_dir+args.file_prefix+".dat"]
    cmd += ["-p", sep_output_dir+args.file_prefix+".ped"]
    cmd += ["-s", sep_output_dir+args.file_prefix+".snps"]
    cmd += ["-h", sep_output_dir+args.file_prefix+".haplos"]
    cmd += ["--prefix", sep_output_dir+args.file_prefix]
    cmd += ["--crossovermap", sep_output_dir+args.file_prefix+".rec"]
    cmd += ["--errormap", sep_output_dir+args.file_prefix+".erate"]    
elif args.mode == 'read':
    print "I'm reading imputed data from chr %s"%(args.chr)
    if cfgs.get('IMPUTE','sample_data')==None or cfgs.get('IMPUTE','sample_data')=='self':
        #I = ImputeData(snpdata_dir+cfgs.get('IMPUTE','reference_data')%args.chr)
        I = ImputeData(output_dir+"filtered_ref_chr%s.csv"%args.chr)
    else:
        I = ImputeData(output_dir+"filtered_samp_chr%s.csv"%args.chr)
    I.Chromosome_divideup(snp_bin_size, snp_bin_size/10, snp_bin_size/5*2)
    I.MACHinput(sep_output_dir+str(args.chr))
    I.imputed_output(output_dir+'imputed_chr%s'%args.chr)
elif args.mode == 'readdose':
    print "I'm reading imputed data from chr %s"%(args.chr)
    if cfgs.get('IMPUTE','sample_data')==None or cfgs.get('IMPUTE','sample_data')=='self':
        #I = ImputeData(snpdata_dir+cfgs.get('IMPUTE','reference_data')%args.chr)
        I = ImputeData(output_dir+"filtered_ref_chr%s.csv"%args.chr)
    else:
        #I = ImputeData(output_dir+"filtered_samp_chr%s.csv"%args.chr)
        J = ImputeData(output_dir+"filtered_ref_chr%s.csv"%args.chr)
    #I.Chromosome_divideup(snp_bin_size, snp_bin_size/10, snp_bin_size/5*2)
    J.Chromosome_divideup(snp_bin_size, snp_bin_size/10, snp_bin_size/5*2)
    #I.MACHdoseinput(sep_output_dir+str(args.chr))
    J.MACHdoseconvert(sep_output_dir+str(args.chr))
    #I.imputed_output(output_dir+'imputed_chr%s.dose'%args.chr)    
    J.imputed_output(output_dir+'ref_chr%s.dose'%args.chr)
elif args.mode == 'error':
    ' Not used now '
    pass
