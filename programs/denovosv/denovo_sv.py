'''
Created on Dec 5, 2011
@author: dazhe.meng

Retrieve the list of contigs/segments that does not map to thaliana and check whether they can map to lyrata
'''

# important imports
import subprocess, sys
import os
sys.path.append("/home/GMI/dazhe.meng/PyLibrary/")
import scaffolds

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("scaff",help="file containing scaffold seqs", default='')
parser.add_argument("blast",help="file containing blast to reference", default='')
parser.add_argument("-o","--outfolder",help="folder to output to", default='./')
parser.add_argument("-f","--outfilename",help="output file name, default to guessing from scaffold file name", default='')
parser.add_argument("--min_score",help="do not consider hits that falls below this threshold", type=int, default=200)
parser.add_argument("--min_homo",help="do not consider hits that are not at least this percent matched", type=int, default=80)
parser.add_argument("--min_bp_score",help="for a bp (or pair of bp), the minimum score increase to consider", type=int, default=400)
parser.add_argument("--max_homo",help="allow some microhomology at both ends, but cut it at some threshold", type=int, default=30)
parser.add_argument("--max_penalty",help="penalty for wrong chromosome/order, don't really need to set", type=int, default=1000000)
parser.add_argument("--min_N",help="when using scaffolds(contain Ns), minimum to ", type=int, default=30)
parser.add_argument("--min_gap",help="smallest gap to try to squeeze in another alignment", type=int, default=400)
parser.add_argument("--min_event_size",help="minimum size of event to report", type=int, default=100)
parser.add_argument("--min_event_size_ends",help="file containing scaffold seqs", type=int, default=400)
parser.add_argument("-d","--debugctg",help="DEBUG",default='')
args = parser.parse_args()

if args.outfilename == '':
    args.outfilename = args.scaff.split('/')[-1].split('.')[0]

S = scaffolds.scaffolds()
print "Reading scaffolds"
S.read_scafseq(args.scaff)
print "Reading blast information"
S.load_blastres(args.blast,'at')
print "Computing break points"
S.get_sv(debugctg=args.debugctg)
S.sort_and_output_sv(args.outfolder+args.outfilename+'.raw')
