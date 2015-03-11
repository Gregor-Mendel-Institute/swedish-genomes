'''
Created on Mar 7, 2012
@author: dazhe.meng

An attempt to parse the new files de novo sv calling

Brainstorming section:
types of big SVs:
plenvar: reference insertion
alseq: alternative sequence (insertion together with deletion) 
sdbp: same direction break point (deletion and/or transloc)
inbp: inverted break point (inversion w/wo transloc)

all directions must be standardized
due to the directional data now available for breakpoints, the 'pos' member no longer always refer to the lower coordinate
instead, now each breakpoint always has 2 pos, with direction fixed --> bp1      bp2 <--
'''

import sys, os

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-d","--ddir",help="directory containing raw sv output", default='/projects/1001genomes/dazhe/denovo_sv/new/')
parser.add_argument("-o","--outfn",help="output file name", default='test.out')
parser.add_argument("--mappingallowance",help="",type=int,default=200)
parser.add_argument("--nallowance",help="",type=int,default=250)
parser.add_argument("--longdistthres",help="",type=int,default=10000)
parser.add_argument("--output_minalfreq",help="",type=int,default=1)
parser.add_argument("--merge_traceback_limit",help="",type=int,default=8)
args = parser.parse_args()

mappingallowance = args.mappingallowance # allow 25 bp overlap between alignments
nallowance = args.nallowance # allow N  
longdistthres = args.longdistthres
output_minalfreq = args.output_minalfreq
merge_traceback_limit = args.merge_traceback_limit

class sv():
    ' A class for describing SV events '
    def __init__(self, chrl, posl, chrr, posr, acc, type, contig_dist, ref_dist, intn):
        self.chrl = chrl
        self.posl = posl # a standard position, defaults to the lower coordinates for breakpoints
        self.chrr = chrr
        self.posr = posr
        self.acc = acc
        self.type = type
        self.cdist = contig_dist
        self.rdist = ref_dist
        self.intn = intn
        
    def __cmp__(self, other):
        if self.chrl == other.chrl:
            return cmp(self.posl, other.posl)
        else:
            return cmp(self.chrl, other.chrl)
    
    def __str__(self):
        return "%s:%s.%s;%s.%s;%s;%s;%s;%s"%(self.type,self.chrl,self.posl,self.chrr,self.posr,self.acc,self.cdist,self.rdist,self.intn)

# Now use a unified dictionary to store all svs
svtypes = ['plenvar','alseq','sdbp','inbp','tbp']
l_sv = {'plenvar':[],'alseq':[],'sdbp':[],'inbp':[],'tbp':[]}



chrdict = {'Chr1':1,'Chr2':2,'Chr3':3,'Chr4':4,'Chr5':5,'chloroplast':6,'mitochondria':7}

print "reading sv..."
listdir = os.listdir(args.ddir)
for fn in listdir:
    if fn.endswith('.raw'):
        f = open(args.ddir+fn)
        a = fn[:-4] # used as accession name
        for l in f:
            o = l.strip().split(',')
            if o[0]=='pplenvar':
                continue # dont parse these yet
            chrl = chrdict[o[1]]
            posl = int(o[2])
            chrr = chrdict[o[3]]
            posr = int(o[4])
            intn = int(o[5])
            cdist = int(o[6])
            if o[7].isdigit():
                rdist = int(o[7])
            else:
                rdist = 'NA'
            l_sv[o[0]].append(sv(chrl, posl, chrr, posr, a, o[0], cdist, rdist, intn))
print "  read %s svs, among these:"%(sum([len(l_sv[t]) for t in svtypes]))
for t in svtypes:
    print "    %s are %s"%(len(l_sv[t]), t)

# sort then group events together
print "sorting..."
for svtype in svtypes:
    l_sv[svtype].sort()

print "grouping similar events:"
l_grouped = []
class sv_group(sv):
    def __init__(self, firstsv):
        self.chrl = firstsv.chrl
        self.posl = firstsv.posl # a standard position, defaults to the lower coordinates for breakpoints
        self.chrr = firstsv.chrr
        self.posr = firstsv.posr
        self.type = firstsv.type
        self.member = [firstsv]
        if firstsv.intn!=0: 
            self.al = nallowance*1.414 # no 'clean' member
            self.fuzzy = True
        else:
            self.al = nallowance
            self.fuzzy = False
    
    def add_member(self, newmem):
        n = len(self.member)
        self.posl = (self.posl*n+newmem.posl)/(n+1)
        self.posr = (self.posr*n+newmem.posr)/(n+1)
        self.member.append(newmem)
  
    def sim(self, osv): # now very strict, except for when there are Ns between the bps
        if self.chrl == osv.chrl and self.chrr == osv.chrr:
            if abs(self.posl-osv.posl)<mappingallowance*2 and abs(self.posr-osv.posr)<mappingallowance*2:
                self.add_member(osv)
                return True
        return False
    
    def roughsim(self, osv): # the 'rough' version for svs, to allow for the cases when there are n's
        if self.chrl == osv.chrl and self.chrr == osv.chrr:
            if abs(self.posl-osv.posl)<self.al and abs(self.posr-osv.posr)<self.al:
                self.add_member(osv)
                return True
        return False    
    
    def check_allele_freq(self):
        ' run this after all members are found '
        self.accs = list(set([m.acc for m in self.member]))
        self.paccs = list(set([m.acc if m.intn==0 else -1 for m in self.member])-set([-1])) # 'pure' or 'clean' accs
        self.accs.sort()
        return (len(self.accs), len(self.member), len(self.paccs))

    def recalculate_dist(self):
        ' run this after all members are found, recalculate the average rdist and cdist '
        if not self.fuzzy:
            l_cdist = []
            for mem in self.member:
                if mem.intn==0:
                    l_cdist.append(mem.cdist)
            try:
                self.cdist = sum(l_cdist)/len(l_cdist)
            except:
                print "Error calculating cdist: %s"
                for mem in self.member:
                    print "  %s"%mem
                self.cdist = 0
        else:
            self.cdist = sum([mem.cdist for mem in self.member])/len(self.member)
        if self.chrl == self.chrr:
            if self.type == 'tbp':
                self.rdist = 'transloc'
            else:
                self.rdist = self.posr-self.posl
        else:
            self.rdist = 'interchr'

for svtype in svtypes:
    print "Merging %s"%svtype
    begcount = len(l_grouped)
    l_traceback = [] # to store svs that have Ns and thus checked later
    print "  round 1: clean svs"
    for sv in l_sv[svtype]:
        if sv.intn != 0:
            l_traceback.append(sv)
        else:
            matched = 0
            for g in l_grouped[len(l_grouped)-merge_traceback_limit:len(l_grouped)]:
                if g.sim(sv):
                    matched += 1
            if not matched:
                l_grouped.append(sv_group(sv))
            elif matched>1:
                sys.stderr.write('Multiple match detected: %s\n'%sv)
    # in the second round, try to fit them with existing calls first, then create new ones
    print "  round 2: %s sv with n"%(len(l_traceback))
    newgroups = []
    min_g_ind = begcount # sets minimum index to improve search efficiency
    min_ng_ind = 0 
    for sv in l_traceback:
        matched = 0
        # clean matches
        if min_g_ind<len(l_grouped):
            while l_grouped[min_g_ind].chrl<sv.chrl or l_grouped[min_g_ind].posl<sv.posl-l_grouped[min_g_ind].al:
                min_g_ind+=1
                if min_g_ind>=len(l_grouped):
                    break
            g_ind = min_g_ind
            while g_ind<len(l_grouped):
                if l_grouped[g_ind].chrl>sv.chrl or l_grouped[g_ind].posl>sv.posl+l_grouped[g_ind].al:
                    break # terminate search
                if l_grouped[g_ind].roughsim(sv):
                    matched +=1
                g_ind +=1
        # rough matches
        if min_ng_ind<len(newgroups):
            while newgroups[min_ng_ind].chrl<sv.chrl or newgroups[min_ng_ind].posl<sv.posl-newgroups[min_ng_ind].al:
                min_ng_ind+=1
                if min_ng_ind>=len(newgroups):
                    break
            ng_ind = min_ng_ind
            while ng_ind<len(newgroups):
                if newgroups[ng_ind].chrl>sv.chrl or newgroups[ng_ind].posl>sv.posl+newgroups[ng_ind].al:
                    break # terminate search
                if newgroups[ng_ind].roughsim(sv):
                    matched +=1
                ng_ind +=1            
        # create new groups if not matched
        if not matched:
            newgroups.append(sv_group(sv))
        #elif matched>1:
        #    sys.stderr.write('Multiple match detected: %s\n'%sv)
    l_grouped+=newgroups
    print "  Done, %s groups found out of %s total"%(len(l_grouped)-begcount, len(l_sv[svtype]))

# output 1
print "Outputting raw SVs"
l_grouped.sort()
f = open(args.ddir+args.outfn,'w')
for sv in l_grouped:
    freq,nfreq,cfreq = sv.check_allele_freq()
    sv.recalculate_dist()
    if freq >= output_minalfreq:
        #oo = sv.check_nonsweep()
        f.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(sv.type,sv.chrl,sv.posl,sv.chrr,sv.posr,cfreq,freq,nfreq,sv.cdist,sv.rdist,','.join(sv.accs)))
