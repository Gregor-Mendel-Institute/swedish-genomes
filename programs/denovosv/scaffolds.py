'''
Created on Dec 5, 2011
@author: dazhe.meng

Retrieve the list of contigs/segments that does not map to thaliana and check whether they can map to lyrata
'''

# important imports
import subprocess, sys
import collections
import os
#sys.path.append("/home/GMI/dazhe.meng/PyLibrary/")
from optparse import OptionParser
from operator import itemgetter
import numpy as np
import bisect
#import globals

scaf = collections.namedtuple('scaf', ['cov','len','seq','blast_at','blast_ly'])
blastentry = collections.namedtuple('blastentry', ['target','qfrom','qto','tfrom','tto','hom','score','strand'])

# sv call paramters

class scaffolds():
    'Class to read and intepret scaffolds'
    def __init__(self):
        self.scafs = {} # dictionary to store scaffold data, use old dict for compatibility
        self.scafl = []  
        self.blast = {} # dictionary to store blast results, each library having a separate sub-dictionary

    def read_scafseq(self, scafseqfile):
        'read soap scafseq files'
        f = open(scafseqfile)
        o = f.readline().strip().split()
        name = o[0][1:]
        cov = float(o[1])
        seq = ""
        for l in f:
            if l[0]=='>':
                self.scafs[name]=scaf(cov, len(seq), seq, [], [])
                self.scafl.append(name)
                seq = ""
                o = l.strip().split()
                name = o[0][1:]
                cov = float(o[1])
            else:
                seq += l.strip()
        self.scafs[name]=scaf(cov, len(seq), seq, [], [])
        self.scafl.append(name)       
    
    def load_blastres(self, blastresfile, whichlib='at'):
        'read blast output'
        f = open(blastresfile)
        for l in f:
            o = l.strip().split()
            hom = float(o[2])
            scafname = o[0]
            target = o[1]
            qfrom = int(o[6])
            qto = int(o[7])
            tfrom = int(o[8])
            tto = int(o[9])
            score = int(float(o[11]))
            strand = 1 if tfrom < tto else -1
            if strand == 1:
                newentry = blastentry(target, qfrom, qto, tfrom, tto, hom, score, strand)
            else:
                newentry = blastentry(target, qfrom, qto, tto, tfrom, hom, score, strand)
            eval('self.scafs[scafname].blast_%s.append(newentry)'%whichlib)

    def check_scaffold(self, whichscaf, thres = 0):
        A_at = self.get_covarray(whichscaf, 'at')
        A_ly = self.get_covarray(whichscaf, 'ly')
        len_nocov = 0
        len_covwly = 0
        for i in xrange(len(A_at)):
            if A_at[i]!=-1: # Not N
                if A_at[i]<=thres:
                    len_nocov +=1
                    if A_ly[i]>thres:
                        len_covwly +=1
        return (len_nocov, len_covwly)

    def get_covarray(self, whichscaf, whichlib):
        ss = self.scafs[whichscaf]
        blasts = eval('ss.blast_%s'%whichlib)
        A = [0]*ss.len
        for i, char in enumerate(ss.seq):
            if char == 'N':
                A[i] = -1
        for hit in blasts:
            for i in xrange(hit.qfrom-1, hit.qto):
                A[i] = max(A[i], hit.hom)
        return A

    def check_all(self, thres = 0):
        tlen_nocov = 0
        tlen_covwly = 0
        count = 0
        for scafname in self.scafl:
            len_nocov, len_covwly = self.check_scaffold(scafname, thres)
            tlen_nocov += len_nocov
            tlen_covwly += len_covwly
            count += 1
            if count % 1000 == 0:
                if tlen_nocov == 0:
                    ratio = float(tlen_covwly/tlen_nocov)*100
                else:
                    ratio = 0.0
                print "processed %s scaffolds, current ratio is at %s/%s (%.2f%%)"%(count, tlen_covwly, tlen_nocov, ratio)        
        return (tlen_nocov, tlen_covwly)

    def get_sv(self, min_score=200, min_homo=80, min_bp_score=400, max_homo = 30, max_penalty = 1000000, min_N = 30, min_gap = 400, out_f='test_sv.csv', min_event_size = 100, min_event_size_ends=400, debugctg=''):
        '''
        calls structural variants from data, based on dynamic programming best collinear alignments on the reference
        
        have to:
        - preferentially put things in the same region together
        
        details: 
        S(i,j) = Best score up to alignment hit i, that ends in j (j<=i)
        path(j) = corresponding path (0,1 array of whether some alignments are taken)
        min_score: do not consider hits that falls below this threshold
        min_bp_score: for a bp (or pair of bp), the minimum score increase to consider  
        max_homo: allow some microhomology at both ends, but cut it at some threshold
                  also should be smaller than min_event_size!
        max_penalty: maximum allowing penalty for a SV
        min_event_size: minimum size of event to report
        '''
        self.break_points = []
        count = 0
        
        for scaf in self.scafs:
            if debugctg and scaf != debugctg:
                continue
            print "now processing: %s"%scaf
            # sort the blast hits in order of increasing qbegin, and initialize the score/path matrix
            l_hit = []
            for hit in self.scafs[scaf].blast_at:
                if hit.score>=min_score and hit.hom>=min_homo:
                    l_hit.append(hit)
            l_hit.sort(key=itemgetter(1))
            n_hit = len(l_hit)
            if n_hit == 0: # no significant hit at all! skip
                continue
            S = [0]*n_hit # have bigger matrix just make coding easier, and hopefully not too much overhead!
            path = np.zeros((n_hit,n_hit),dtype=np.int8)
            len_scaf = self.scafs[scaf].len
            if len_scaf > max_penalty:
                max_penalty = len_scaf

            for i, hit in enumerate(l_hit):
                print 'hit%s:\t%s-%s\t%s:%s-%s\t%s'%(i, hit.qfrom, hit.qto, hit.target, hit.tfrom, hit.tto, hit.strand)

            # obtain the position of Ns in the scaffold, with some minimum length
            l_Nsegb = []
            l_Nsege = []
            cur_pos = 0
            seq = self.scafs[scaf].seq + 'E' # last E is for the case when the last letter is N, which should not happen!
            Nseed = 'N'*min_N
            Nind = seq[cur_pos:].find(Nseed)
            while Nind!=-1:
                i = Nind+cur_pos
                while seq[i]=='N':
                    i+=1
                l_Nsegb.append(Nind+cur_pos)
                l_Nsege.append(i)
                cur_pos = i
                Nind = seq[cur_pos:].find(Nseed)
            
            def count_n(pfrom, pend):
                pfrom -= 1
                if pend < pfrom:
                    return 0
                ibeg = bisect.bisect_right(l_Nsegb, pfrom)
                iend = bisect.bisect_left(l_Nsege, pend)
                sum = 0
                for i in xrange(ibeg,iend):
                    sum += l_Nsege[i]-l_Nsegb[i]
                return sum
                
            # dynamic programming/greedy algorithm to choose the best linear alignment:
            #    1: with extreme penalty for crossing huge genomic distance, find the best linear alignment
            #    2: now allowing for discordant joins one (with its pair) at a time, go though the other alignments within 'gaps' 
            for i in xrange(n_hit):
                l_score = []
                for j in xrange(i):
                    if l_hit[j].qto - l_hit[i].qfrom > max_homo: # if the overlap is too great
                        neighbor_score = -max_penalty # just some arbitrary huge penalty
                    else:
                        if l_hit[i].target != l_hit[j].target:
                            neighbor_score = -max_penalty
                        elif l_hit[i].strand != l_hit[j].strand:
                            neighbor_score = -max_penalty
                        else: 
                            if l_hit[i].strand == 1 and l_hit[i].tfrom < l_hit[j].tto:
                                neighbor_score = -max_penalty
                            elif l_hit[i].strand == -1 and l_hit[i].tto > l_hit[j].tfrom:
                                neighbor_score = -max_penalty
                            else:
                                neighbor_score = max(min(-(l_hit[i].tfrom-l_hit[j].tto),-(l_hit[i].qfrom-l_hit[j].qto)),-max_penalty)
                    l_score.append(S[j]+neighbor_score)
                l_score.append(-max_penalty)
                max_j = np.argmax(l_score)
                S[i] = l_score[max_j]
                if S[i]< -(l_hit[i].qfrom-1):
                    S[i] = -(l_hit[i].qfrom-1) # meaning taking i as the first hit
                else:
                    path[i] = path[max_j]
                S[i]+=l_hit[i].score
                path[i,i] = 1
                print l_score
                print S
                print path[i]
                print
            for i in xrange(n_hit):
                S[i] -= len_scaf - l_hit[i].qto # meaning taking i as the last hit
            max_i = np.argmax(S)
            self.cur_path = []
            for i in xrange(n_hit):
                if path[max_i][i]==1:
                    self.cur_path.append(i)
                    
            # iterative steps, might look mad!
            def close_gap(i_hitb=None, i_hite=None, i_cur_path=0):
                print "closing gap: %s-%s"%(i_hitb, i_hite)
                if i_hitb == None: # closing begining
                    p_beg = 1
                    p_end = l_hit[i_hite].qfrom
                    i_hitb = -1
                elif i_hite == None:
                    p_end = len_scaf
                    p_beg = l_hit[i_hitb].qto
                    i_hite = n_hit
                else:
                    if i_hite-i_hitb==1:
                        return 1 #means closed
                    p_beg = l_hit[i_hitb].qto
                    p_end = l_hit[i_hite].qfrom
                    if p_end-p_beg < min_gap:
                        return 1
               
                print "range: %s-%s"%(p_beg, p_end)

                i_hitb +=1 # python range begin from i_hitb
                while i_hitb<n_hit: 
                    if p_beg-l_hit[i_hitb].qfrom>max_homo:
                        i_hitb+=1
                    else:
                        break
                nn_hit = i_hite - i_hitb
                if nn_hit<=0:
                    return 1 
                S = [0]*nn_hit 
                path = np.zeros((nn_hit,nn_hit),dtype=np.int8)
                for i in xrange(nn_hit):
                    l_score = []
                    for j in xrange(i):
                        if l_hit[j+i_hitb].qto - l_hit[i+i_hitb].qfrom > max_homo: # if the overlap is too great
                            neighbor_score = -max_penalty # just some arbitrary huge penalty
                        else:
                            if l_hit[i+i_hitb].target != l_hit[j+i_hitb].target:
                                neighbor_score = -max_penalty
                            elif l_hit[i+i_hitb].strand != l_hit[j+i_hitb].strand:
                                neighbor_score = -max_penalty
                            elif l_hit[i+i_hitb].tfrom < l_hit[j+i_hitb].tto:
                                neighbor_score = -max_penalty
                            else:
                                neighbor_score = max(min(-(l_hit[i+i_hitb].tfrom-l_hit[j+i_hitb].tto),-(l_hit[i+i_hitb].qfrom-l_hit[j+i_hitb].qto)),-max_penalty)
                        l_score.append(S[j]+neighbor_score)
                    l_score.append(-max_penalty)
                    max_j = np.argmax(l_score)
                    S[i] = l_score[max_j]
                    if S[i]< -(l_hit[i+i_hitb].qfrom-p_beg):
                        S[i] = -(l_hit[i+i_hitb].qfrom-p_beg) # meaning taking i as the first hit
                    else:
                        path[i] = path[j]
                    S[i]+= l_hit[i+i_hitb].score
                    path[i,i] = 1
                for i in xrange(nn_hit):
                    S[i] -= max(0,p_end - l_hit[i+i_hitb].qto) # meaning taking i as the last hit
                    if l_hit[i+i_hitb].qto - p_end > max_homo:
                        S[i] -= max_penalty # same rule as above
                if max(S)+(p_end-p_beg)<min_bp_score: # improvement not significant
                    return 1
                max_i = np.argmax(S)
                sub_path = []
                for i in xrange(nn_hit):
                    if path[max_i][i]==1:
                        sub_path.append(i+i_hitb)   
                self.cur_path = self.cur_path[:i_cur_path]+sub_path+self.cur_path[i_cur_path:] 
                return 0
            
            # close the beginning
            test_beg = close_gap(i_hite=self.cur_path[0], i_cur_path=0)
            while not test_beg:
                test_beg = close_gap(i_hite=self.cur_path[0], i_cur_path=0)
                
            # close the end
            test_end = close_gap(i_hitb=self.cur_path[-1],i_cur_path=len(self.cur_path))
            while not test_end:
                test_end = close_gap(i_hitb=self.cur_path[-1],i_cur_path=len(self.cur_path))            
            
            # close the middle
            i_close_gap = 0
            while i_close_gap != len(self.cur_path)-1:
                test_close = close_gap(i_hitb=self.cur_path[i_close_gap],i_hite=self.cur_path[i_close_gap+1],i_cur_path=i_close_gap+1)
                while not test_close:
                    test_close = close_gap(i_hitb=self.cur_path[i_close_gap],i_hite=self.cur_path[i_close_gap+1],i_cur_path=i_close_gap+1)
                i_close_gap+=1
                    
            # finally! have the final sequence, now output the variations
            # OUTPUT FORMAT: type,lchr,lpos,rchr,rpos,N,qdist,tdist
            # types: pplenvar: partial plenvar
            #        plenvar:
            #        inbp
            #        sdbp
            #        alseq
            # lchr, lpos, rchr, rpos: positions of the breakpoints, smallest first. Note that, for pplenvar and plenvar, only one set is used!
            # N: number of Ns in the sequence between the breakpoints
            # size: size of the 
            
            # check beginning for possible insertion
            cur_path=self.cur_path

            # sanity check
            print "total len: %s"%len_scaf
            for i in cur_path:
                hit = l_hit[i]
                print 'hit%s:\t%s-%s\t%s:%s-%s\t%s'%(i, hit.qfrom, hit.qto, hit.target, hit.tfrom, hit.tto, hit.strand)

            hit = l_hit[cur_path[0]]
            if hit.qfrom>min_event_size_ends+1:
                if hit.strand == 1:
                    self.break_points.append(('pplenvar',hit.target,hit.tfrom,'-','-','NA',-hit.qfrom+1,0,scaf)) # the sign indicate direction
                else:
                    self.break_points.append(('pplenvar',hit.target,hit.tto,'-','-','NA',hit.qfrom-1,0,scaf))
            # check ending for possible insertion
            hit = l_hit[cur_path[-1]]
            if len_scaf-hit.qto>min_event_size_ends:
                if hit.strand == 1:
                    self.break_points.append(('pplenvar',hit.target,hit.tto,'-','-','NA',len_scaf-hit.qto,0,scaf)) # the sign indicate direction
                else:
                    self.break_points.append(('pplenvar',hit.target,hit.tfrom,'-','-','NA',-(len_scaf-hit.qto),0,scaf))
            # check the middle for any event
            for i in xrange(len(cur_path)-1):
                hit1 = l_hit[cur_path[i]]
                hit2 = l_hit[cur_path[i+1]]
                if hit1.strand != hit2.strand: # inverted break point
                    '''
                    if hit1.target>hit2.target or (hit1.target==hit2.target and hit1.tfrom>hit2.tfrom,scaf):
                        hitt2 = hit1
                        hitt1 = hit2
                    else:
                        hitt1 = hit1
                        hitt2 = hit2
                    tdist = ['interchr','interchr']
                    if hit1.target == hit2.target:
                        tdist = [hitt2.tto-hitt1.tto, hitt2.tfrom-hitt1.tfrom]
                    if hit1.strand == 1:
                        self.break_points.append(('inbp',hitt1.target,hitt1.tto,hitt2.target,hitt2.tto,count_n(hit1.qto, hit2.qfrom),max(hit2.qfrom-hit1.qto,0),tdist[0],scaf))
                    else:
                        self.break_points.append(('inbp',hitt1.target,hitt1.tfrom,hitt2.target,hitt2.tfrom,count_n(hit1.qto, hit2.qfrom),max(hit2.qfrom-hit1.qto,0),tdist[1],scaf))
                    '''
                    if hit1.strand == 1: # the +- case, sequences are to the left of breakpoints -> smaller coord first
                        obp1 = (hit1.target, hit1.tto)
                        obp2 = (hit2.target, hit2.tto)
                        ohits = [obp1,obp2]
                        ohits.sort()
                        if ohits[0][0]!=ohits[1][0]:
                            tdist = 'interchr'
                        else:
                            tdist = ohits[1][1]-ohits[0][1]
                        self.break_points.append(('inbp',ohits[0][0],ohits[0][1],ohits[1][0],ohits[1][1],count_n(hit1.qto, hit2.qfrom),hit2.qfrom-hit1.qto,tdist,scaf))
                    else: # the -+ case, sequences are to the right of breakpoints -> larger coord first
                        obp1 = (hit1.target, hit1.tfrom)
                        obp2 = (hit2.target, hit2.tfrom)
                        ohits = [obp1,obp2]
                        ohits.sort()
                        if ohits[0][0]!=ohits[1][0]:
                            tdist = 'interchr'
                        else:
                            tdist = ohits[0][1]-ohits[1][1]
                        self.break_points.append(('inbp',ohits[1][0],ohits[1][1],ohits[0][0],ohits[0][1],count_n(hit1.qto, hit2.qfrom),hit2.qfrom-hit1.qto,tdist,scaf))
                else:
                    if hit1.target != hit2.target:
                        if hit1.target>hit2.target:
                            hitt2 = hit1
                            hitt1 = hit2
                        else:
                            hitt1 = hit1
                            hitt2 = hit2
                        if hit1.strand == 1:
                            self.break_points.append(('sdbp',hitt1.target,hitt1.tto,hitt2.target,hitt2.tfrom,count_n(hit1.qto, hit2.qfrom),max(hit2.qfrom-hit1.qto,0),'interchr',scaf))
                        else:
                            self.break_points.append(('sdbp',hitt1.target,hitt1.tfrom,hitt2.target,hitt2.tto,count_n(hit1.qto, hit2.qfrom),max(hit2.qfrom-hit1.qto,0),'interchr',scaf))
                    else:
                        qdist = max(hit2.qfrom - hit1.qto,0)
                        ncount = count_n(hit1.qto, hit2.qfrom)
                        target = hit1.target # the same thing, so doesn't matter
                        if hit1.strand == 1:
                            pos1 = hit1.tto
                            pos2 = hit2.tfrom
                        else:
                            pos1 = hit2.tto
                            pos2 = hit1.tfrom
                        if pos1-pos2>max_homo:
                            self.break_points.append(('tbp',target,pos1,target,pos2,ncount,qdist,'transloc',scaf))
                        tdist = max(pos2-pos1,0)
                        if tdist>=min_event_size or qdist>=min_event_size:
                            if tdist<max_homo: # insertion
                                self.break_points.append(('plenvar',target,pos1,target,pos2,ncount,qdist,tdist,scaf))
                            elif qdist<max_homo: # deletion/sd transloc
                                self.break_points.append(('sdbp',target,pos1,target,pos2,ncount,qdist,tdist,scaf))
                            else:
                                if float(ncount)/qdist > 0.8: # set a hard thres now, should change to param later!
                                    continue
                                self.break_points.append(('alseq',target,pos1,target,pos2,ncount,qdist,tdist,scaf))
            count +=1
            if count % 100 == 0:
                print "processed %s scaffolds"%count
                
    def sort_and_output_sv(self, outfile='test_sv.csv', flag_big=False, flag_slocus=False):
        print "sorting..."
        self.break_points.sort(key=itemgetter(1,2))
        print "outputting to %s..."%outfile
        f = open(outfile,'w')
        if flag_big:
            f2 = open(outfile+'.big','w') # testing only
        if flag_slocus:
            f3 = open(outfile+'.ac','w') # for s locus a c breakpoint
            
        for sv in self.break_points:
            f.write(','.join(map(str,sv))+'\n')
            if flag_big:
                if sv[0] == 'tbp' or sv[0] == 'inbp' or sv[7] == 'interchr':
                    f2.write(','.join(map(str,sv))+'\n')
                else:
                    if sv[7]>500000:
                        f2.write(','.join(map(str,sv))+'\n')
            if flag_slocus:
                if sv[1]!=sv[3]:
                    if sv[1]!='-' and sv[3]!='-':
                        f3.write(','.join(map(str,sv))+'\n') 

    def get_novel_contig(self, outputfile, min_len=0):
        'identifies novel contigs (and those that are not even in ly)'
        novels = []
        for scafname in self.scafl:
            ss = self.scafs[scafname]
            if ss.len<min_len:
                continue
            if len(ss.blast_at)==0:
                novels.append(scafname)
        fo = open(outputfile,'w')
        folist = open(outputfile+'.list','w')
        for novelscaf in novels:
            S = self.scafs[novelscaf]
            fo.write('>%s\n'%novelscaf)
            folist.write('%s\t%s\n'%(novelscaf,S.len))
            for i_chunk in xrange(S.len/100+1): # divide into smaller chunks
                fo.write(S.seq[i_chunk*100:min(S.len,(i_chunk+1)*100)])
            fo.write('\n')

    def get_partial_novel_contig(self, outputfile, min_len=300, thres=80):
        'identifies and extract contigs that contain a significant novel part'
        novels = []
        for scafname in self.scafl:
            ss = self.scafs[scafname]
            if len(ss.blast_at)==0:
                novels.append(scafname)
            else:
                at_no_match = 0
                A_at = self.get_covarray(scafname, 'at')
                for i in xrange(len(A_at)):
                    if A_at[i]!=-1: # Not N
                        if A_at[i]<=thres:
                            at_no_match+=1
                if at_no_match >= min_len:
                    novels.append(scafname)
        fo = open(outputfile,'w')
        folist = open(outputfile+'.list','w')
        for novelscaf in novels:
            S = self.scafs[novelscaf]
            fo.write('>%s\n'%novelscaf)
            folist.write('%s\t%s\n'%(novelscaf,S.len))
            for i_chunk in xrange(S.len/100+1): # divide into smaller chunks
                fo.write(S.seq[i_chunk*100:min(S.len,(i_chunk+1)*100)])
            fo.write('\n')

    def get_all_novel_contig(self, outputfile, min_len=300, thres=80):
        'identifies and extract significant novel part of contigs (and do not output the parts that are linked to at!)'
        novels = []
        fo = open(outputfile,'w')
        folist = open(outputfile+'.list','w')
        for scafname in self.scafl:
            ss = self.scafs[scafname]
            if len(ss.blast_at)==0:
                if ss.len>=min_len:
                    fo.write('>%s\n'%scafname)
                    folist.write('%s\t%s\tNA\n'%(scafname,ss.len))
                    for i_chunk in xrange(ss.len/100+1): # divide into smaller chunks
                        fo.write(ss.seq[i_chunk*100:min(ss.len,(i_chunk+1)*100)])
                    fo.write('\n')
            else:
                A_at = self.get_covarray(scafname, 'at')
                at_loc = None
                output_index = 0
                last_no_match = None
                for i in xrange(len(A_at)):
                    if A_at[i]!=-1 and A_at[i]<=thres:
                        if last_no_match == None:
                            last_no_match = i    
                    else:
                        if last_no_match != None:
                            if i-last_no_match > min_len: # strict greater because i is actually the index of the first base that fails
                                # an eligible piece found! now output
                                if not at_loc:
                                    at_loc = (ss.blast_at[0].target, ss.blast_at[0].tfrom)
                                fo.write('>%s_%s\n'%(scafname,output_index))
                                outputseq = ss.seq[last_no_match:i]
                                outlen = len(outputseq)
                                folist.write('%s_%s\t%s\t%s:%s\n'%(scafname,output_index,outlen,at_loc[0],at_loc[1]))
                                for i_chunk in xrange(outlen/100+1): # divide into smaller chunks
                                    fo.write(outputseq[i_chunk*100:min(outlen,(i_chunk+1)*100)])
                                fo.write('\n')
                                output_index+=1
                            last_no_match = None

