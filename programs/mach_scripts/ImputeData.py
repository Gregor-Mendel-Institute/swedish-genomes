"""
This class script convert data from what I have to the forms as required by specific programmes
It can mask different subsets of the data to behave like it's hidden

Despite the numbering, this comes after ImputeData3 and is intended to be a final version of the ImputeData class

Now, since vardata is much more cleaned up and efficient, a lot of code in this class has much more streamlined than before

Major Updates:
05242011 1. Updating the accession hiding part to current expectation. The old set_hidden_accs and related stuff are removed. 
         In this version, the hiding is achieved by modifying the data matrix itself while storing a copy of hidden part, which I believe is a mucher cleaner way then modifying individual output functions
         2. Adding Output and Input functions to all the various program that Im going to try
"""

# Some temporary stuff which should be moved into class variables later
MACH_samp_snps = "test_machsamp.snps"
MACH_haplo_snps = "test_machref180.snps"
MACH_crossover_map = "test_mach.rec"
MACH_error_map = "test_mach.erate"
BEAGLE_markers = "test_ref180g.bmarkers"

import vardata
from copy import deepcopy
import random
import subprocess
import numpy as np
# optional imports, might move to individual parts
import gzip

class ImputeData:
    def __init__(self, genotype_file, program = "NPUTE", window_size = 100, check_accuracy = False, missing_allele='?'):
        self.genotype = vardata.vardata(genotype_file,loadformat=1) #loadformat was 0
        self.sampgeno = None # For programs that demand separate genotype file
        self.program = program
        self.imputed_genotype = None
        self.window_size = window_size
        self.checking_accuracy = check_accuracy
        self.hidden_acc_data = []
        self.hidden_acc_ind = None
        self.sample_snp_inds = [] # This stores the indices of snps in sample data
        self.num_imputed = 0

    # Input functions, some are reimplemented from the old ImputeData3 class (and probably should die)
    def add_imputing_accs(self, imputing_snp_file, start_ind=0, end_ind=0, replacing = True, inserting = False):
        ' Add the accs to be imputed (hence, normally inserting should be false)' 
        acc_merge_subset = 'all' if start_ind==0 and end_ind==0 else '%s:%s'%(start_ind, end_ind)
        merge_stat = self.genotype.merge_from_file(imputing_snp_file, add_snps = inserting, add_accs = True, replace_allele = replacing, acc_merge_subset = acc_merge_subset)
        print "added %s accessions"%(merge_stat[1])
        print "added %s SNPs and changed status of %s SNPs"%(merge_stat[0], merge_stat[4])
        self.num_imputed = merge_stat[1]

    def set_imputing_snps(self, imputing_snp_file):
        ' Set part of the snps to be considered the imputing set (This is just infill with replacing = False and inserting = False)'
        merge_stat = self.genotype.merge_from_file(imputing_snp_file, add_snps = False, add_accs = False, replace_allele = False)
        print "%s markers among %s are now considered imputing"%(merge_stat[4], self.genotype.num_var)
        if merge_stat[3]!=0:
            print "There is %s/%s allele inconsistencies, or in fraction %s"%(merge_stat[2], merge_stat[3], float(merge_stat[2])/merge_stat[3])
        else:
            print "No overlapping allele information"
    
    def infill_imputing_snps(self, imputing_snp_file, replacing = True, inserting = True):
        ' Insert snps from the imputing file if they dont already exist (This is just add_imputing_accs without start_ind and end_ind)'
        merge_stat = self.genotype.merge_from_file(imputing_snp_file, add_snps = inserting, add_accs = True, replace_allele = replacing)
        print "added %s SNPs and changed status of %s SNPs"%(merge_stat[0], merge_stat[4])
        if merge_stat[3]!=0:
            print "There is %s/%s allele inconsistencies, or in fraction %s"%(merge_stat[2], merge_stat[3], float(merge_stat[2])/merge_stat[3])
        else:
            print "No overlapping allele information"

    def filter_snps(self, **kwargs):
        ' Filter the SNP data, or rather, just calls the filter_snps function in vardata '
        self.genotype.filter_snps(**kwargs)
        
    def separate_data(self):
        ' Separate the imputing and sample data into two '
        self.target_genotype = vardata.vardata()
        self.target_genotype.accessions = self.genotype.accessions[-self.num_imputed:]
        self.target_genotype.vars = self.genotype.vars
        self.target_genotype.data_matrix = self.genotype.data_matrix[:,-self.num_imputed:]
        self.target_genotype.data_format = self.genotype.data_format
        try:
            self.genotype.remove_acc((len(self.genotype.accessions)-self.num_imputed,len(self.genotype.accessions)))
        except:
            print len(self.genotype.accessions)
            print self.num_imputed
            print self.genotype.data_matrix.shape

    # Cross-validation related functions
    def __get_sample_snp_inds(self):
        if not self.sample_snp_inds:
            for i, var in enumerate(self.genotype.vars):
                if var.indicator==1:
                    self.sample_snp_inds.append(i)
 
    def hide_acc(self, accind, remove_data=True, missing_in_samp_geno=False):
        ' this updated function hide one accessions in the data (and store them in separate data) '
        if self.hidden_acc_data: # Note: this checks the data, so if data is restored, the index can still be used to find which acc was hidden
            print "Cannot hide accession %s: Already hiding accession %s!"%(accind, self.hidden_acc_ind)
            return
        if accind<0 or accind>=self.genotype.num_acc:
            raise ValueError("Cannot hide accession %s: Index out of bound"%accind)
        else:
            self.__get_sample_snp_inds()
            self.hidden_acc_ind = accind
            self.hidden_acc_data = list(self.genotype.data_matrix[:,accind])
            if remove_data: # for programs that requires separate files for input_output
                # Construct separate genotype data, notice that I should set up a copy functionin vardata
                self.sampgeno = vardata.vardata() # initialize a new class
                self.sampgeno.accessions = [self.genotype.accessions[accind]]
                if not missing_in_samp_geno:
                    self.sampgeno.data_matrix = np.empty((len(self.sample_snp_inds), 1),dtype=self.genotype.matrix_dtype)
                    for i, snp_i in enumerate(self.sample_snp_inds):
                        self.sampgeno.vars.append(deepcopy(self.genotype.vars[snp_i]))
                        self.sampgeno.vars[-1].ind = i
                        self.sampgeno.data_matrix[i,0] = self.genotype.data_matrix[self.genotype.vars[snp_i].ind,accind]
                    self.sampgeno.num_var = len(self.sample_snp_inds)
                else:
                    self.sampgeno.vars = deepcopy(self.genotype.vars)
                    self.sampgeno.data_matrix = np.tile(np.array(self.genotype.miss_allele,dtype=self.genotype.matrix_dtype),(len(self.genotype.data_matrix), 1))
                    for i in self.sample_snp_inds:
                        ii = self.genotype.vars[i].ind
                        self.genotype.data_matrix[ii,0] = self.hidden_acc_data[ii] 
                    self.sampgeno.num_var = self.genotype.num_var                   
                self.sampgeno.num_acc = 1
                # Remove corresponding accession from genotype  
                self.genotype.remove_acc([accind])              
            else:
                self.genotype.data_matrix[:,accind]=[self.genotype.miss_allele]*len(self.genotype.data_matrix)
                for i in self.sample_snp_inds:
                    ii = self.genotype.vars[i].ind
                    self.genotype.data_matrix[ii,accind] = self.hidden_acc_data[ii]
            
    def restore_hidden_acc_data(self):
        self.genotype.data_matrix[:,self.hidden_acc_ind] = self.hidden_acc_data
        self.hidden_acc_data = []
        
    def cross_validate_check(self):
        ' Cross validation check, currently for only 1 accessions and does not consider missing alleles '
        self.error_check = []
        for i in xrange(len(self.hidden_acc_data)):
            if self.hidden_acc_data[i] == self.imputed_hidden_acc_data[i]:
                self.error_check.append(0)
            else:
                self.error_check.append(1)
        
    def accuracy_table(self, outputname):
        ' checks for accuracy of imputation and output it in a csv '
        error_file = open(outputname, "w")
        error_file.write("Chromosome,Position,"+",".join(self.hidden_accessions)+"\n")
        error_count = [0]*(self.genotype.acclen+1) # Keep the count of the errors, probably dont need so many but nevermind
        indl = self.genotype.index_acc(self.hidden_accessions)
        total_checked = 0
        for i in xrange(len(self.genotype.data)):
            true_entry = self.genotype.data[i]
            if true_entry.indicator == 1: # Do not test for the imputing snps
                continue
            impute_entry = self.imputed_genotype.data[i]
            if true_entry != impute_entry:
                print true_entry
                print impute_entry
                raise ValueError("Imputed and target entry positions not matching!")
            error_file.write("%s,%s"%(true_entry.chr, true_entry.pos))
            num_error = 0
            for i in indl:
                if true_entry.data[i]=='?' or true_entry.data[i]=='M':
                    error_file.write(",?")
                elif true_entry.data[i]!=impute_entry.data[i]:
                    num_error += 1
                    total_checked += 1
                    error_file.write(",1")
                else:
                    total_checked += 1
                    error_file.write(",0")
            error_file.write("\n")
            error_count[num_error]+=1     
        #print error_count
        total_error = sum([i*error_count[i] for i in range(len(error_count))])
        print float(total_error)/total_checked
        return error_count      
    
    # Some functions to set additional parameters
    def set_temp_file(self, tmpfileprefix = ""):
        ' Set the name of the temporary files if needed '
        self.tmpfileprefix = tmpfileprefix
        
    def set_program_dir(self, program_dir):
        self.program_dir = program_dir

    def set_info_file_folder(self, info_folder = "/home/GMI/nordborg-group-common/dazhe/Test_SNPdata/"):
        self.info_file_folder = info_folder
        
    def set_additional_params(self, **kwargs):
        self.additional_params = kwargs

    # Wrapper interface functions with imputation programs    
    def output(self):
        ' Wrapper function that calls the slave functions '
        if self.program == "NPUTE":
            self.NPUTEoutput(self.tmpfileprefix+"_IN") 
        elif self.program == "BEAGLE":
            self.BEAGLEoutput(self.tmpfileprefix)
        elif self.program == "MACH":
            self.MACHoutput(self.tmpfileprefix)
 
    def input(self):
        ' Wrapper function that calls the slave functions '
        if self.program == "NPUTE":
            self.NPUTEinput(self.tmpfileprefix+"_OUT")
        elif self.program == "BEAGLE":
            self.BEAGLEinput(self.tmpfileprefix)
        elif self.program == "MACH":
            self.MACHinput(self.tmpfileprefix)

    def run(self):
        ' Functions that calls the imputing program and returns the retcode '
        if self.program == "NPUTE":
            inputf_name = self.tmpfileprefix+"_IN"
            outputf_name = self.tmpfileprefix+"_OUT"
            for chunkind in xrange(len(self.chunks)):
                retcode = subprocess.call("%sNPUTE.py -w %s -i %s_%s -o %s_%s -f 2"%(self.program_dir, self.window_size, inputf_name, chunkind, outputf_name, chunkind), shell=True)
                #retcode = subprocess.call("%sNPUTE.py -w %s -i %s_%s -o %s_%s"%(self.program_dir, self.window_size, inputf_name, chunkind, outputf_name, chunkind), shell=True)
            return retcode
        elif self.program == "BEAGLE":
            cmd = ["java","-Xmx4000m","-jar", "%sbeagle.jar"%self.program_dir]
            cmd.append("phased=%s"%(self.tmpfileprefix+".ref"))
            cmd.append("phased=%s"%(self.tmpfileprefix+".samp"))
            cmd.append("markers=%s"%(self.info_file_folder+BEAGLE_markers))
            cmd.append("missing=?")
            cmd.append("omitprefix=true")
            cmd.append("out=%s"%(self.tmpfileprefix))
            retcode = subprocess.call(cmd)
        elif self.program == "MACH":
            if self.additional_params['mode'] == "run":
                cmd = ["%smach1"%self.program_dir, "--greedy", "--rounds %s"%self.additional_params['rounds']]
                cmd += ["-d", self.tmpfileprefix+".dat"]
                cmd += ["-p", self.tmpfileprefix+".ped"]
                cmd += ["-s", self.tmpfileprefix+".snps"]
                cmd += ["-h", self.tmpfileprefix+".haplos"]
                cmd += ["--prefix", self.tmpfileprefix]
                retcode = subprocess.call(cmd)
            elif self.additional_params['mode'] == "run_noref":
                cmd = ["%smach1"%self.program_dir, "--rounds %s"%self.additional_params['rounds']]
                cmd += ["-d", self.tmpfileprefix+".dat"]
                cmd += ["-p", self.tmpfileprefix+".ped"]
                cmd += ["--prefix", self.tmpfileprefix]
                retcode = subprocess.call(cmd)                
            elif self.additional_params['mode'] == "mle":
                cmd = ["%smach1"%self.program_dir, "--greedy", "--mle"]
                cmd += ["-d", self.tmpfileprefix+".dat"]
                cmd += ["-p", self.tmpfileprefix+".ped"]
                cmd += ["-s", self.tmpfileprefix+".snps"]
                cmd += ["-h", self.tmpfileprefix+".haplos"]
                cmd += ["--prefix", self.tmpfileprefix]
                cmd += ["--crossovermap", self.tmpfileprefix+".rec"]
                cmd += ["--errormap", self.tmpfileprefix+".erate"]                
    
    def cleanup(self):
        if self.program == "NPUTE":
            self.NPUTEcleanup()
    
    # Actual interface functions with specific programs
    def NPUTEoutput(self, filename):
        ' This is really easy, as NPUTE format is almost the same '
        for chunkind in xrange(len(self.chunks)): # works for datasets with all chromosomes only
            print "Outputting chunk %s/%s"%(chunkind+1, len(self.chunks))
            outfile = open(filename+"_%s"%chunkind,"w")
            for entry in self.genotype.vars[self.chunks[chunkind][0]:self.chunks[chunkind][1]]:
                tmpentrydata = self.genotype.data_matrix[entry.ind]     
                outfile.write(",".join(tmpentrydata)+"\n")
    
    def NPUTEinput(self, filename):
        if self.checking_accuracy == True:
            #self.imputed_genotype = deepcopy(self.genotype) # This eats memory, maybe just use it for comparison
            self.imputed_hidden_acc_data = [self.genotype.miss_allele]*len(self.hidden_acc_data)
            read_index = 0
            last_chunk_end = 0
            for chunkind in xrange(len(self.chunks)):
                print "Reading chunk %s/%s"%(chunkind+1, len(self.chunks))
                infile = open(filename+"_%s"%chunkind)
                if last_chunk_end > self.chunks[chunkind][0]: # There is overlap
                    read_index -= self.chunkbuffer
                    for bufferind in xrange(self.chunkbuffer):
                        infile.readline()
                last_chunk_end = self.chunks[chunkind][1]
                for line in infile:
                    if line.strip()=="":
                        break
                    self.imputed_hidden_acc_data[read_index] = line.strip().split(",")[self.hidden_acc_ind]
                    read_index+=1
            infile.close()
        else:
            self.imputed_genotype = self.genotype
            read_index = 0
            last_chunk_end = 0
            for chunkind in xrange(len(self.chunks)):
                print "Reading chunk %s/%s"%(chunkind+1, len(self.chunks))
                infile = open(filename+"_%s"%chunkind)
                if last_chunk_end > self.chunks[chunkind][0]: # There is overlap
                    read_index -= self.chunkbuffer
                    for bufferind in xrange(self.chunkbuffer):
                        infile.readline()
                last_chunk_end = self.chunks[chunkind][1]
                for line in infile:
                    if line.strip()=="":
                        break
                    self.imputed_genotype.data_matrix[read_index] = line.strip().split(",")
                    read_index+=1
            infile.close()
        #subprocess.call("rm -rf "+filename+"_%s"%chunkind, shell=True) # should implement better clean-up function

    def NPUTEcleanup(self):
        for chunkind in xrange(len(self.chunks)):
            inputf_name = self.tmpfileprefix+"_IN_%s"
            outputf_name = self.tmpfileprefix+"_OUT_%s"
            subprocess.call("rm -f "+inputf_name+"_%s"%chunkind, shell=True)
            subprocess.call("rm -f "+outputf_name+"_%s"%chunkind, shell=True)

    # Following is mostly temporary code for cross-validation
    def MACHoutput(self, filename, mode = 'Single'):
        if mode == 'Trial':
            self.genotype.output(filename+".haplos", format="MACH_haplo")
            self.sampgeno.output(filename+".ped", format="Merlin")
        elif mode == 'Single':
            for chunkind in xrange(len(self.chunks)): # works for datasets with all chromosomes only
                print "Outputting chunk %s/%s"%(chunkind+1, len(self.chunks))
                self.genotype.output(filename+"_%s"%chunkind+".dat", "MACH_markers", startind=self.chunks[chunkind][0] , endind=self.chunks[chunkind][1])
                self.genotype.output(filename+"_%s"%chunkind+".ped", "Merlin", startind=self.chunks[chunkind][0] , endind=self.chunks[chunkind][1])         
        elif mode == 'Double':
            for chunkind in xrange(len(self.chunks)):
                print "Outputting chunk %s/%s"%(chunkind+1, len(self.chunks))
                self.genotype.output(filename+"_%s"%chunkind+".snps", "MACH_haplomarkers", startind=self.chunks[chunkind][0] , endind=self.chunks[chunkind][1])
                self.genotype.output(filename+"_%s"%chunkind+".haplos", "MACH_haplo", startind=self.chunks[chunkind][0] , endind=self.chunks[chunkind][1])         
                self.target_genotype.output(filename+"_%s"%chunkind+".dat", "MACH_markers", startind=self.chunks[chunkind][0] , endind=self.chunks[chunkind][1])
                self.target_genotype.output(filename+"_%s"%chunkind+".ped", "Merlin", startind=self.chunks[chunkind][0] , endind=self.chunks[chunkind][1])         
    
    def MACHinput(self, filename): # Needs some serious testing
        if self.checking_accuracy == True:
            self.imputed_hidden_acc_data = [self.genotype.miss_allele]*len(self.hidden_acc_data)
            f_dose = gzip.open(filename+".mldose.gz")
            f_allele = open(filename+".mlinfo")
            doses = f_dose.readline().strip().split()[2:]
            f_allele.readline()
            for i, dosestr in enumerate(doses):
                dosage = float(dosestr)
                alleles = f_allele.readline().strip().split()[1:3]
                self.imputed_hidden_acc_data[i] = alleles[0] if dosage>=1 else alleles[1]
        else:
            self.imputed_genotype = self.genotype
            read_index_start = 0
            last_chunk_end = 0
            for chunkind in xrange(len(self.chunks)):
                print "Reading chunk %s/%s"%(chunkind+1, len(self.chunks))
                f_dose = gzip.open(filename+"_%s"%chunkind+".dose.gz")
                f_allele = open(filename+"_%s"%chunkind+".info")
                # get allele information first
                f_allele.readline()
                l_allele = []
                for line in f_allele:
                    l_allele.append(line.strip().split()[1:3])
                # now get genotype info line by line
                if last_chunk_end > self.chunks[chunkind][0]: # There is overlap
                    read_index_start -= self.chunkbuffer
                    offset = self.chunkbuffer
                else:
                    offset = 0
                last_chunk_end = self.chunks[chunkind][1]
                for accind, line in enumerate(f_dose):
                    doses = line.strip().split()[2:]
                    read_index = read_index_start
                    for snpind, dosage in enumerate(doses[offset:]):
                        al = l_allele[snpind+offset][0] if float(dosage) >= 1 else l_allele[snpind+offset][1]
                        self.imputed_genotype.data_matrix[read_index, accind]= al
                        read_index+=1
                read_index_start = read_index 

    def MACHdoseinput(self, filename): # try dosage file for Matt
        self.imputed_genotype = self.genotype
        self.imputed_genotype.data_matrix = np.zeros(shape=self.imputed_genotype.data_matrix.shape, dtype=np.dtype('a5'))
        read_index_start = 0
        last_chunk_end = 0
        for chunkind in xrange(len(self.chunks)):
            print "Reading chunk %s/%s"%(chunkind+1, len(self.chunks))
            f_dose = gzip.open(filename+"_%s"%chunkind+".dose.gz")
            f_allele = open(filename+"_%s"%chunkind+".info")
            # get allele information first
            f_allele.readline()
            l_allele = []
            for line in f_allele:
                l_allele.append(line.strip().split()[1:3])
            # now get genotype info line by line
            if last_chunk_end > self.chunks[chunkind][0]: # There is overlap
                read_index_start -= self.chunkbuffer
                offset = self.chunkbuffer
            else:
                offset = 0
            last_chunk_end = self.chunks[chunkind][1]
            for accind, line in enumerate(f_dose):
                doses = line.strip().split()[2:]
                read_index = read_index_start
                for snpind, dosage in enumerate(doses[offset:]):
                    #al = l_allele[snpind+offset][0] if float(dosage) >= 1 else l_allele[snpind+offset][1]
                    self.imputed_genotype.data_matrix[read_index, accind]= '%.3f'%float(dosage)
                    read_index+=1
            read_index_start = read_index

    def MACHdoseconvert(self, filename): # try dosage file for Matt
        self.imputed_genotype = self.genotype
        tmp_data_matrix = np.zeros(shape=self.imputed_genotype.data_matrix.shape, dtype=np.dtype('a5'))
        read_index_start = 0
        last_chunk_end = 0
        for chunkind in xrange(len(self.chunks)):
            print "Reading chunk %s/%s"%(chunkind+1, len(self.chunks))
            f_allele = open(filename+"_%s"%chunkind+".info")
            # get allele information first
            f_allele.readline()
            l_allele = []
            for line in f_allele:
                l_allele.append(line.strip().split()[1:3])
            # now get genotype info line by line
            if last_chunk_end > self.chunks[chunkind][0]: # There is overlap
                read_index_start -= self.chunkbuffer
                offset = self.chunkbuffer
            else:
                offset = 0
            last_chunk_end = self.chunks[chunkind][1]
            read_index = read_index_start
            for snpid, lal in enumerate(l_allele[offset:]):
                #al = l_allele[snpind+offset][0] if float(dosage) >= 1 else l_allele[snpind+offset][1]
                for accind in xrange(self.imputed_genotype.num_acc): # a weird inefficient way of doing this?
                    orii = self.imputed_genotype.data_matrix[read_index, accind]
                    tmp_data_matrix[read_index, accind]= '0.000' if orii == lal[1] else '2.000'
                read_index+=1
            read_index_start = read_index
        self.imputed_genotype.data_matrix = tmp_data_matrix

    def MACHcleanup(self):
        pass
    
    def BEAGLEoutput(self, filename):
        self.genotype.output(filename+".ref", format="BEAGLE")
        self.sampgeno.output(filename+".samp", format="BEAGLE")
    
    def BEAGLEinput(self, filename):
        if self.checking_accuracy == True:
            self.imputed_hidden_acc_data = [self.genotype.miss_allele]*len(self.hidden_acc_data)
            f_imputed = gzip.open(filename+".samp.phased.gz")
            f_imputed.readline()
            for i, l in enumerate(f_imputed):
                self.imputed_hidden_acc_data[i] = l.strip().split()[2]
    
    def BEAGLEcleanup(self):
        pass

    def fastPHASEoutput(self, filename):
        self.genotype.output(filename+".haplo", format="fastPHASE_haplo")
        self.sampgeno.output(filename+".diplo", format="fastPHASE_diplo")
    
    def fastPHASEinput(self, filename):
        if self.checking_accuracy == True:
            self.imputed_hidden_acc_data = [self.genotype.miss_allele]*len(self.hidden_acc_data)
            f_imputed = open(filename+".samp.phased.gz")
            f_imputed.readline()
            for i, l in enumerate(f_imputed):
                self.imputed_hidden_acc_data[i] = l.strip().split()[2]
    
    def fastPHASEcleanup(self):
        pass

    # Very useful utility functions
    def Chromosome_divideup(self, chunksize = 10000, chunkbuffer = 500, minchunksize = 2000):
        ' Divides up the genome file into small chunks; the genome file should contain single chromosomes only '
        self.chunks = []
        self.chunkbuffer = chunkbuffer
        if len(self.genotype.vars)<chunksize+chunkbuffer:
            self.chunks.append([0, len(self.genotype.vars)])
        else:
            num_chunk = len(self.genotype.vars)/chunksize
            self.chunks.append([0, chunksize+chunkbuffer])
            for i in xrange(1, num_chunk):
                self.chunks.append([chunksize*i-chunkbuffer, chunksize*(i+1)+chunkbuffer])
            if len(self.genotype.vars)-num_chunk*chunksize<minchunksize:
                self.chunks[-1][1] = len(self.genotype.vars) # Extend to the end if the remaining chunk is too small
            else:
                self.chunks.append([chunksize*num_chunk - chunkbuffer, len(self.genotype.vars)]) # use the remainder as the last chunk
    
    def Output_chunks(self, filename):
        ' Output the number of chunks Chromosome_divideup'
        fo = open(filename,'w')
        fo.write('%s'%len(self.chunks))
    
    # Output functions
    def imputed_output(self, outfilename, include_imputing_snp = True, include_imputing_acc = False):
        ' Output the imputed accessions '
        f = open(outfilename,"w")
        if include_imputing_acc == False:
            index_start = -self.num_imputed
        else:
            index_start = 0
        f.write("Chromosome,Position,"+",".join(self.genotype.accessions[index_start:])+"\n")            
        for var in self.imputed_genotype.vars:
            if var.indicator == 1 and include_imputing_snp == False:
                continue
            else:
                f.write("%s,%s,"%(var.chr, var.pos)+",".join(self.imputed_genotype.data_matrix[var.ind][index_start:])+"\n")

    def error_output(self, outfilename, include_imputing_snp=False):
        ' Output cross validation error file, currently only working for 1 acc '
        f = open(outfilename,"w")
        for var in self.genotype.vars:
            if var.indicator == 1 and include_imputing_snp == False:
                continue
            else:
                f.write(str(self.error_check[var.ind])+"\n")

    # Debugging functions
    def Simulate_Error(self, rate = 0.05):
        self.imputed_genotype = deepcopy(self.genotype)
        for snp in self.imputed_genotype.data:
            snp.data = ["X" if s!="?" and random.random()<rate else s for s in snp.data]
