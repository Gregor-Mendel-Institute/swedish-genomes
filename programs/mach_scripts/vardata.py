"""
This is a completely new data structure, and despite the number, is actually offspring of vardata3
The goal is to improve vardata3's terrible memory efficiency by using Numpy's int8 matrix
varentry is kept somewhat as a header information for retrieving real data

Note: It might be interesting to link the main class in each varentry

Major Updates:

"""

import tools # For bufcount
import time
from bisect import bisect_left, bisect_right
from copy import deepcopy
import numpy as np

class varentry:
    def __init__(self, chr, pos, data_ind = -1, moreposinfo = None, indicator = 0):
        self.ind = data_ind #data index for main class
        self.chr = chr
        self.pos = pos
        self.completepos = moreposinfo
        self.indicator = indicator
        
    def __eq__(self, other):
        return self.chr == other.chr and self.pos == other.pos
    
    def __cmp__(self, other):
        ' A direct comparison function to compare pos1itions, should be later revamped to individual functions for faster processing '
        if self.chr == other.chr:
            return cmp(self.pos, other.pos)
        else:
            return cmp(self.chr, other.chr)
    
    def __sub__(self, other):
        ' subtraction function, for finding position difference '
        if self.chr == other.chr:
            return self.pos - other.pos
        else:
            return None

class vardata:
    """
    This class is (supposedly) designed to encompass all types of variation data: SNP or indels (CNV)
    accessions (x): accession list
    data (x*y): [chr, pos]+[list of chars]
    2010.03.01 overhauled readfromfile with inclusion of ability to compound genotypes (for plink tped) and display progress
    2010.05.24 reworking this to really fit all types of variation dataset, with the basic architecture more or less intact; this heavily involves varentry, a new class to retain position and genotype information. It will also now use heap exclusively for fast retrieval and storage
    2010.05.27 now uses bisect instead of heapq
    2011.05.12 reworked, uses numpy int8 table now
    """
    def __init__(self, filename="", loadformat = 0, dataformat=0, min_maf = 0, max_misrate = 1, max_allele = 2):
        self.accessions = []
        self.vars = []
        self.dict = {'A':'A','T':'T','G':'G','C':'C','-':'-','NA':'?','AG':'?','CT':'?','AT':'?','AC':'?','CG':'?','GT':'?', '?':'?', 'R':'?', 'Y':'?', 'S':'?', 'W':'?', 'K':'?', 'M':'?', 'B':'?', 'D':'?', 'H':'?', 'V':'?', 'N':'?'}
        self.revdict = {'A':'A','T':'T','G':'G','C':'C','-':'-','?':'NA', 'R':'NA', 'Y':'NA', 'S':'NA', 'W':'NA', 'K':'NA', 'M':'NA', 'B':'NA', 'D':'NA', 'H':'NA', 'V':'NA', 'N':'NA'}
        self.num_acc = 0
        self.max_allele = 2
        self.load_format = loadformat
        self.data_format = dataformat
        self.bool_reordered = False # stores if the data has been changed in order
        self.G = None # For storing genome information
        self.allele_info = None
        if filename!="":
            self.readfromfile(filename)
 
    def readfromfile(self, filename, displayprogress=1, sep=",", compound_genotype=0, header=1, miss_allele='?'):
        'Initialize from a file, 1 for bjarni format'
        fi = open(filename)
        if displayprogress==1:
            print "Reading variation data from "+filename
            starting_time = time.time()
        if header == 1:
            line = fi.readline()
            u = line.strip().split(sep)
            self.accessions = u[2:]
            self.num_acc = len(self.accessions)
        self.read_sep = sep
        self.read_compound_genotype = compound_genotype
        
        # Initializing data
        total_num_lines = tools.bufcount(filename)
        if self.load_format in [0,1]:
            self.num_var = total_num_lines - header
        if self.data_format == 0:
            self.miss_allele = miss_allele
            self.matrix_dtype = np.dtype('a1')
        elif self.data_format == 1:
            self.miss_allele = -1
            self.matrix_dtype = np.int8
        self.data_matrix = np.empty((self.num_var,self.num_acc),self.matrix_dtype) # string array
        self.__read_index = 0 # stores the index for the next entry of data
        
        if displayprogress == 1:
            print "%s used for line counting and memory allocation"%(time.time()-starting_time)
            pd = tools.progress_display(total_num_lines)
            if header:
                pd.inc()
            for line in fi:
                new_var, new_data = self.__parse_data_entry__(line)
                self.data_matrix[self.__read_index] = new_data
                self.vars.append(new_var)
                self.__read_index+=1
                pd.inc()
        else:
            for line in fi:
                new_var, new_data = self.__parse_data_entry__(line)
                self.data_matrix[self.__read_index] = new_data
                self.vars.append(new_var)
                self.__read_index+=1
        if displayprogress == 1:
            print "Loading complete, %s entries loaded"%len(self.vars) 
        # sanity check
        if self.num_var!=len(self.vars):
            raise ValueError("Should have %s, but only %s loaded!"%(self.num_var, len(self.var))) 

    def merge_from_file(self, filename, add_snps = False, add_accs = True, replace_allele = False, indicator = 1, sep=",", compound_genotype=0, header=1, acc_merge_subset = 'all', use_ref_for_missing = False, reorder = "True"):
        """
        merges snpdata from another snpdata file
        
        return format: (NUM SNPS ADDED, NUM ACCS ADDED, NUM SNP ALLELES REPLACED/IN DISAGREEMENT, NUM SNP ALLELES OBSERVED IN BOTH, NUM OF ORIGINAL SNPS THAT HAS ITS INDICATOR CHANGED)
        
        parameters:
        - add_snps: whether new snps from the file should be appended
        - replace_allele: whether to replace alleles that are in disagreement with the new file; default only replacing those that are missing in the original set
        - add_accs: whether to append new accessions to snp array
        - indicator: change the indicator of snps that are in the new file
        - acc_merge_subset: specify a subset of the accessions to be merged (this was designed for imputation) in the format of startind:endind 
        - use_ref_for_missing: specify whether to use the reference allele to fill in the alleles that were not included in one of the data sets 
        - reorder: whether to automatically sort the var list after merging
        
        potential problems:
        - SNP data is stored as an array, and thus when expanding the data it is necessary to allocate the same amount of memory again
        - A search is performed for each entry in the merging file. While i tried to optimize it by skipping over all previously searched space, it might still increase computation time greatly when the merging file is large
        
        assumptions:
        - original data is already sorted (should be)
        - new data file is sorted and in Bjarni's format
        """
        if not add_snps and not add_accs and not replace_allele and indicator == 0:
            raise ValueError("Something is wrong, you dont seem to want to change anything in the old data!")
        else:
            print "Merging data..."
        
        fi = open(filename)
        if header == 1:
            line = fi.readline()
            u = line.strip().split(sep)
            new_accessions = u[2:]
            self.read_sep = sep
        else:
            raise IOError("New snp file does not contain accession information!")
        self.read_compound_genotype = compound_genotype
        
        # if appending accessions, expand the data accordingly
        if add_accs == True:
            if acc_merge_subset != 'all':
                acc_subset_start_ind, acc_subset_end_ind = [int(tmp_i) for tmp_i in acc_merge_subset.split(':')]
                subset_new_accs = new_accessions[acc_subset_start_ind:acc_subset_end_ind]
            else:
                subset_new_accs = new_accessions
            num_added_acc = self.expand_acc(subset_new_accs, use_ref_for_missing)
            print "Added new accessions"
        else:
            num_added_acc = 0

        # index the accessions in the new snp data set
        valid_indices = []
        new_acc_indices = self.index_acc(new_accessions)
        for i in xrange(len(new_acc_indices)):
            if new_acc_indices[i]!= -1:
                valid_indices.append(i)
        
        # initialize lists to store new snp information, they are simply unused if add_snps = False
        l_new_vars = []
        l_new_data = []
        
        # initialize variables to store dataset difference information
        n_examined = 0
        n_errors = 0
        n_changed_indicator = 0
        self.retrieve_init()
        
        # the main iterative part of this function
        print "Reading merging data file..."
        for line in fi:
            new_var, new_data = self.__parse_data_entry__(line)
            existing = self.retrieve_next(new_var)
            if existing: # this snp was in the old set
                existing.indicator = indicator
                n_changed_indicator += 1
                for overlapped_ind in valid_indices:
                    existing_allele = self.data_matrix[existing.ind][new_acc_indices[overlapped_ind]]
                    if existing_allele == "?":
                        self.data_matrix[existing.ind][new_acc_indices[overlapped_ind]] = new_data[overlapped_ind]
                    else:
                        n_examined += 1
                        if existing_allele != new_data[overlapped_ind]:
                            n_errors += 1
                            if replace_allele == True:
                                self.data_matrix[existing.ind][new_acc_indices[overlapped_ind]] = new_data[overlapped_ind]
            else:
                if add_snps:
                    new_var.indicator = indicator
                    l_new_vars.append(new_var)
                    if use_ref_for_missing:
                        tmp_data_row = np.tile(np.array(self.G.get_pos(new_var.chr, new_var.pos), dtype=self.matrix_dtype), self.num_acc)
                    else:
                        tmp_data_row = np.tile(np.array(self.miss_allele,dtype=self.matrix_dtype),(self.num_acc))
                    for overlapped_ind in valid_indices:
                        tmp_data_row[new_acc_indices[overlapped_ind]] = new_data[overlapped_ind]
                    l_new_data.append(tmp_data_row)
                    self.__read_index += 1
        
        # finally, adding snps if they are good
        if add_snps:
            if len(l_new_vars)==0:
                print "No SNP to add"
            else:
                print "Adding SNPs"
                self.data_matrix = np.append(self.data_matrix, l_new_data, axis=0)
                self.vars += l_new_vars
                self.num_var += len(l_new_vars)
                if reorder:
                    self.vars.sort()
                self.bool_reordered = True
            
        return (len(l_new_vars), num_added_acc, n_errors, n_examined, n_changed_indicator)

    def __parse_data_entry__(self, line):
        'Internal slave function to parse the data line for readfromfile'
        u = line.strip().split(self.read_sep)
        completepos = None
        if self.load_format == 1:
            if u[1].isdigit()==True: # A snp
                pos = int(u[1])
            else:
                pos = int(u[1].split(":")[0]) # Correspond to indel position format
                completepos = u[1]
            new_var=varentry(int(u[0]),pos,self.__read_index, completepos)
            new_data=[self.dict[i] for i in u[2:]]
        elif self.load_format == 3: #plink tped
            if self.read_compound_genotype == 0:
                pass # to be fixed
            else: # this is for the case where plink tped file have multiple column for each allele in multiploidal dataset
                new_var=[int(u[0]),int(u[3]),[u[1]]+self.__compound_genotype__(u[4:])]
        elif self.load_format == 4: # 0, 1
            dict01 = {'0':0, '1':1, '?':-1}
            if u[1].isdigit()==True: # A snp
                pos = int(u[1])
            else:
                pos = int(u[1].split(":")[0]) # Correspond to indel position format
                completepos = u[1]
            new_var=varentry(int(u[0]),pos,self.__read_index,completepos) 
            new_data=[dict01[i] for i in u[2:]]
        else:
            if u[1].isdigit()==True: # A snp
                pos = int(u[1])
            else:
                pos = int(u[1].split(":")[0]) # Correspond to indel position format
                completepos = u[1]
            new_var=varentry(int(u[0]),pos,self.__read_index,completepos)   
            new_data=[i for i in u[2:]]
        return (new_var, new_data)

    def __compound_genotype__(self, genotype_list):
        'Internal, compounds a genotype in the case of plink tped format'
        new_genotypelist = []
        if len(genotype_list)%2 != 0:
            raise "Incorrect number of snp data in compound genotype list"
        for i in xrange(len(genotype_list)/2):
            new_genotypelist.append("".join(sorted([genotype_list[i],genotype_list[i+1]])))
        return new_genotypelist
            
    def expand_acc(self, moreaccs, use_ref_as_missing = False):
        ' Expand the accession list to include new accessions from the moreaccs list'
        new_acc_no = len(set(moreaccs)-set(self.accessions))
        if use_ref_as_missing:
            if self.data_format == 1:
                raise ValueError("Add reference allele as genotype currently not working for 0,1 format")
            self.load_genome()
        if new_acc_no == 0:
            print "Nothing new in the accessions"
        else:
            if not use_ref_as_missing:
                matrix_added_part = np.tile(np.array(self.miss_allele,dtype=self.matrix_dtype),(len(self.data_matrix), new_acc_no))
            else:
                matrix_added_part = np.empty((len(self.data_matrix), new_acc_no), self.matrix_dtype)
                for var in self.vars:
                    matrix_added_part[var.ind] = np.repeat(self.G.get_pos(var.chr, var.pos), new_acc_no)
            self.data_matrix=np.append(self.data_matrix,matrix_added_part, axis=1) # This is really a potential place where memory requirement explode! But for now i shall leave it be
            self.num_acc += new_acc_no
            self.accessions += list(set(moreaccs)-set(self.accessions))
        return new_acc_no

    def remove_acc(self, accs, keep_listed = False):
        ' Remove specific accessions from data, here maintaining order of data is very important! '
        new_accessions = []
        
        # Smart detection of argument format
        if isinstance(accs, list):
            acc_index_list = accs
        elif isinstance(accs, tuple):
            acc_index_list = range(accs[0], accs[1])
        elif isinstance(accs, int):
            acc_index_list = [accs]
        elif isinstance(accs, str):
            acc_index_list = [self.index_acc([accs])[0]]
            
        acc_index_list.sort()
        if acc_index_list[-1]>=self.num_acc or acc_index_list[0]<0:
            print "Invalid accession indices: %s"%acc_index_list
            return
        if keep_listed:
            for accind in self.acc_index_list:
                new_accessions.append(self.accessions[accind])
        else:
            for accind in xrange(self.num_acc):
                if accind not in acc_index_list:
                    new_accessions.append(self.accessions[accind])
        self.accessions = new_accessions
        self.num_acc = len(new_accessions)
        if keep_listed:
            new_data_matrix = self.data_matrix.take(acc_index_list, axis=1)
        else:
            new_data_matrix = np.delete(self.data_matrix,acc_index_list, axis=1)
        self.data_matrix = new_data_matrix

    def validate_self(self):
        ' some basic sanity check '
        if self.__read_index != len(self.data_matrix):
            print "Read index unsynced!"
        else:
            print "Read index = %s (sane)"%self.__read_index
        
        if self.num_acc != len(self.accessions):
            print "Number of accessions unsynced!"
        else:
            print "Number of accessions = %s (sane)"%self.num_acc
            
        if self.num_var != len(self.vars):
            print "Number of markers unsynced!"
        else:
            print "Number of markers = %s (sane)"%self.num_var

    def cleanup_deleted(self):
        ' reserved as a function to clean up the data matrix and remove all rows that are unlinked to a snp class'
        pass

    def index_acc(self, accs):
        'Index the entries in accs in self.accessions, missing ones would not be indexed and receive -1!'
        indexlist = []
        for i in accs:
            try:
                indexlist.append(self.accessions.index(i))
            except ValueError:
                indexlist.append(-1)
        return indexlist
    
    def retrieve(self, entry, start_ind = 0):
        ' Test if the dataset already contain an entry with same position '
        try:
            oldentry = self.vars[bisect_left(self.vars[start_ind:], entry)]
        except IndexError:
            return None
        if entry == oldentry:
            return oldentry
        else:
            return None   
    
    def retrieve_init(self):
        self.__retrieve_index = 0
    
    def retrieve_next(self, entry):
        ' Smart retrieve that assumes that querys came sequentially'
        if self.__retrieve_index == self.num_var:
            return None
        else:
            for i in xrange(self.__retrieve_index, self.num_var):
                if self.vars[i]<entry:
                    continue
                elif self.vars[i]==entry:
                    self.__retrieve_index = i
                    return self.vars[i]
                else:
                    self.__retrieve_index = i
                    return None
            self.__retrieve_index = self.num_var
            return None
    
    def iterate(self):
        'iterate through the variation data'
        for i in self.data:
            yield i
    
    def output(self, filename = "", format = 'NoConversion', chr = None, startpos = None, endpos = None, startind = None, endind = None):
        'output the vardata in a format indicated'
        if filename == "":
            filename = "VarData.csv"
        #print "Writing data to %s..."%filename
        var_startind = 0
        var_endind = len(self.data_matrix)
        if chr: # regional output
            var_startind, var_endind = self.get_region(chr, startpos, endpos)
            var2o = self.vars[var_startind:var_endind]
        elif startind != None and endind != None:
            var_startind, var_endind = (startind, endind)
            var2o = self.vars[startind:endind]
        else:
            var2o = self.vars
        
        if format == "Plink": # Transposed Plink format
            tpedfile = open(filename+".tped","w")
            tfamfile = open(filename+".tfam","w")
            for snp in var2o:
                tpedfile.write("%s %s_%s 0 %s "%(snp.chr, snp.chr, snp.pos, snp.pos)+" ".join(["%s %s"%(x, x) for x in self.data_matrix[snp.ind]])+"\n")
            for acc in self.accessions:
                tfamfile.write("%s %s 0 0 1 1\n"%(acc, acc))
        elif format == "fastPHASE_haplo":
            haplofile = open(filename,"w")
            haplofile.write("%s\n"%self.num_acc)
            for accind, accname in enumerate(self.accessions):
                haplofile.write("%s\n"%accname)
                haplofile.write("".join(list(self.data_matrix[:,accind]))+"\n")
        elif format == "fastPHASE_diplo":
            diplofile = open(filename,"w")
            diplofile.write("%s\n"%self.num_acc)
            diplofile.write("%s\n"%self.num_var)
            for accind, accname in enumerate(self.accessions):
                diplofile.write("%s\n"%accname)
                line = "".join(list(self.data_matrix[:,accind]))+"\n"
                diplofile.write(line+line)
        elif format == "Eminim":
            snpdatafile = open(filename+"_snps.tsv","w")
            # accfile = open(filename+"_accs.csv") # not useful when trying
            for entry in var2o:
                snpdatafile.write("%s_%s\t%s\t%s\t%s\t%s\t+\t"%(entry.chr, entry.pos, entry.chr, entry.pos, self.allele_info[entry.ind][0], self.allele_info[entry.ind][1]))
                snpdatafile.write("\t".join([str(al+1) for al in self.data_matrix[entry.ind]]))
                snpdatafile.write("\n")
        elif format == "Merlin": # paired data, is it a problem?
            outfile = open(filename, "w")
            dict_machhaplo = {'A':'A','C':'C','G':'G','T':'T','?':'N','N':'N'}
            for accind, accname in enumerate(self.accessions):
                outfile.write("%s\t%s\t0\t0\tM\t"%(accname, accname))
                outfile.write("\t".join(["{0} {0}".format(dict_machhaplo[al]) for al in self.data_matrix[var_startind:var_endind,accind]])+"\n")            
        elif format == "MACH_haplo": # the haplotype format for MACH, first two columns are ignored, and note ACGT = 1,2,3,4
            outfile = open(filename, "w")
            dict_machhaplo = {'A':'A','C':'C','G':'G','T':'T','?':'N'}  # this file does not allow for missing data!
            for accind, accname in enumerate(self.accessions):
                outfile.write("%s %s "%(accname, accname))
                #outfile.write("".join([dict_machhaplo[al] for al in self.data_matrix[:,accind]])+"\n")
                outfile.write("".join(self.data_matrix[var_startind:var_endind,accind])+"\n")
        elif format == "MACH_markers": # The markers file required by MACH
            outfile = open(filename, "w")
            for entry in var2o:
                outfile.write("M %s_%s\n"%(entry.chr, entry.pos))
        elif format == "MACH_haplomarkers": # The haplotype markers file requried by MACH (yeah, so many different ones)
            outfile = open(filename, "w")
            for entry in var2o:
                outfile.write("%s_%s\n"%(entry.chr, entry.pos))
        elif format == "BEAGLE": # similar to bjarni's but each snp is identified by a single string
            outfile = open(filename, "w")
            outfile.write("I ID "+" ".join(self.accessions)+"\n")
            for entry in var2o:
                outfile.write("M "+" ".join([str(entry.chr)+"_"+str(entry.pos)]+[str(x) for x in self.data_matrix[entry.ind]])+"\n") 
        elif format == "BEAGLEmarkers":
            outfile = open(filename,"w")
            for entry in var2o:
                outfile.write(str(entry.chr)+"_"+str(entry.pos)+"\t%s\t"%entry.pos+"\t".join(list(set(self.data_matrix[entry.ind])-set([self.miss_allele])))+"\n") # Ideally I should have a separate step to extract the alleles            
        elif format == 'Bjarni':
            outfile = open(filename, "w")
            outfile.write("Chromosome, Position,"+",".join(self.accessions)+"\n")
            for entry in var2o:
                outfile.write(",".join([str(entry.chr),str(entry.pos)]+[self.revdict[str(x)] for x in self.data_matrix[entry.ind]])+"\n")
        elif format == 'NoConversion':
            outfile = open(filename, "w")
            outfile.write("Chromosome,Position,"+",".join(self.accessions)+"\n")
            if self.bool_reordered == False and chr == None: # this output only works when not specifying region!
                for ind, entry in enumerate(self.vars):
                    outfile.write(",".join([str(entry.chr),str(entry.pos)]+[str(x) for x in self.data_matrix[ind]])+"\n")  
            else:
                for entry in var2o:
                    outfile.write(",".join([str(entry.chr),str(entry.pos)]+[str(x) for x in self.data_matrix[entry.ind]])+"\n")
    
    def output_sep_chr(self, filename = ""):
        ' output the files in separate chromosomes '
        if filename == "":
            filename = "Vardata_Chr%s.csv"
        current_chr = 0
        for entry in self.vars:
            if entry.chr > current_chr:
                current_chr = entry.chr
                f = open(filename%current_chr, "w")
                f.write("Chromosome,Position,"+",".join(self.accessions)+"\n")
            f.write(",".join([str(entry.chr),str(entry.pos)]+[str(x) for x in self.data_matrix[entry.ind]])+"\n")
        f.close()
 
    def count_var_in_chromosome(self):
        ' count the number of polymorphism sites in each chromosome '
        counts = []
        for entry in self.vars:
            if entry.chr>len(counts):
                counts.append(0)
            counts[entry.chr-1]+=1
        return counts

    def load_genome(self):
        if not self.G:
            import genome
            self.G = genome.Genome()

    # Data manipulation tools
    def convert_into_01(self, store_allele_info=True, force=False):
        if store_allele_info:
            #self.load_genome()
            self.allele_info = np.empty((len(self.data_matrix), 2), dtype=np.dtype("a1"))
        new_data_matrix = np.empty(np.shape(self.data_matrix), dtype = np.int8)
        for varind in xrange(len(self.data_matrix)):
            tmp_l = self.data_matrix[varind]
            l_al = list(set(tmp_l)-set(["?"]))
            if len(l_al)!=2:
                if force==False:
                    print "something wrong at %s: %s alleles!"%(self.vars[varind].pos, len(l_al))
                    print tmp_l
            self.allele_info[varind] = list(set(tmp_l)-set(["?"]))
            d_mini = {"?":-1, self.allele_info[varind][0]:0, self.allele_info[varind][1]:1}
            new_data_matrix[varind] = [d_mini[al] for al in tmp_l]
        self.data_format = 1
        self.data_matrix = new_data_matrix
    
    def get_region(self, chr, startpos, endpos):
        ' return the index list of vars in a certain region; of course, the vars should be in order here '
        start = varentry(chr, startpos)
        end = varentry(chr, endpos)
        left = bisect_left(self.vars, start)
        right = bisect_right(self.vars, end)
        if left == self.num_var or right == 0:
            return ()
        elif left>=right:
            return ()
        else:
            return (left,right)
    
    def filter_snps(self, min_maf = 0, max_mis = 1, max_allele = 0, remove_data = True):
        ' Filters the SNP data according to preset rules, note that SNPs are not actually removed from the data table is there has been reordering of the data_matrix '
        num_filtered = 0
        new_vars = []

        if min_maf == 0 and max_mis == 1 and max_allele  == 0:
            print "Nothing to filter, exiting"
            return
        if self.bool_reordered:
            print "Data table not in order, actual data will not be removed"
        
        if max_allele == 0:
            max_allele = self.num_acc # just so it is the largest possible number of alleles
        for var in self.vars:
            snpdata = list(self.data_matrix[var.ind])
            alleles = list(set(snpdata)-set([self.miss_allele]))
            if len(alleles)>max_allele:
                num_filtered += 1
                continue
            if max_mis <1:
                mis_rate = float(snpdata.count(self.miss_allele))/self.num_acc
                if mis_rate > max_mis:
                    num_filtered += 1
                    continue
            if min_maf >0:
                if len(alleles)>=2:
                    counts = [snpdata.count(i) for i in alleles]
                    counts.sort()
                    maf = float(counts[-2])/sum(counts)
                    if maf < min_maf:
                        num_filtered += 1
                        continue
                else: # no minor allele, consider as 0
                    num_filtered += 1
                    continue
            # every check passed
            new_vars.append(var)
        print "Filtering complete, %s out of %s markers didn't make it"%(num_filtered, self.num_var)
        if num_filtered != 0:
            self.vars = new_vars
            self.num_var = self.num_var - num_filtered
            if remove_data == False:
                self.bool_reordered = True
            else:
                if self.bool_reordered == False: # Note that if it was true, it will remain true
                    print "Removing SNP data from data matrix"
                    keep_indices = []
                    for var_ind, var in enumerate(self.vars): # here it's updated already!
                        keep_indices.append(var.ind)
                        var.ind = var_ind # update the indices to the new matrix that is soon to come
                    self.data_matrix = self.data_matrix.take(keep_indices, axis=0)
                    self.__read_index = len(self.data_matrix)
                    
    def compare_quality(self, other):
        ' Compares the quality of one snp set against another, assuming both are sorted and ordered '
        num_same = 0
        num_diff = 0
        num_missed = 0
        self.retrieve_init()
        l_other_acc_ind = self.index_acc(other.accessions)
        for var in other.vars:
            self_var = self.retrieve_next(var)
            if self_var:
                dat_self = self.data_matrix[self_var.ind]
                dat_other = other.data_matrix[var.ind]
                for acc_ind in xrange(other.num_acc):
                    self_acc_ind = l_other_acc_ind[acc_ind]
                    if self_acc_ind!=-1:
                        if dat_self[self_acc_ind]==self.miss_allele:
                            if dat_other[acc_ind]!=other.miss_allele:
                                num_missed +=1
                        elif dat_self[self_acc_ind]==dat_other[acc_ind]:
                            num_same += 1
                        else:
                            if dat_other[acc_ind]!=other.miss_allele:
                                num_diff += 1
        return (num_same, num_diff, num_missed)
    
    # Some additional functions that serves to calculate LD and such, assumes 01 data format and no missing/maf0 data!
    def countal1(self):
        if self.data_format != 1:
            raise ValueError("SNP data not in 0,1 format!")
        #self.al1counts = np.zeros((len(self.data_matrix)),dtype=np.int32)
        #for i in xrange(len(self.data_matrix)):
        #    self.al1counts[i]=sum(self.data_matrix[i])
        self.al1counts = np.dot(self.data_matrix, np.ones((self.num_acc), dtype=np.int64))
            
    def init_ld(self):
        ' Allocate memory for all tmp variables, and do all necessary conversions '
        if self.data_format == 0:
            print 'Converting data into binary format'
            self.convert_into_01()
        print 'Calculating frequencies'
        self.countal1()
        #print 'Calculating pairwise dot product'
        #self.init_xymat()
        self.x = 0
        self.y = 0
        self.n = self.num_acc
        self.xy = 0
        self.n2D = 0.0
    
    def init_xymat(self):
        self.xym = np.dot(self.data_matrix, np.transpose(self.data_matrix))
    
    def ld(self, ind1, ind2):
        ' calculates r2 ld between two rows of datamatrix with indices ind1 and ind2 '
        self.x = self.al1counts[ind1]
        self.y = self.al1counts[ind2]
        self.xy = sum(self.data_matrix[ind1]*self.data_matrix[ind2])
        #self.xy = self.xym[ind1,ind2]
        self.n2D = float(self.n * self.xy - self.x*self.y)
        return (self.n2D)*(self.n2D)/self.x/self.y/(self.n-self.x)/(self.n-self.y)
