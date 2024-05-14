#!/usr/bin/env python3

import sys
from Bio import  SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import numpy as np
from dataclasses import dataclass
from rep_utils import Kmer_Dict

read_fasta_file = sys.argv[1]  # repeat region read fasta file
blast_outfmt6_file = sys.argv[2] # output of reads (query) blastn against repeat region sequence (subject)
FLANK_LEN = int(sys.argv[3]) # length of flanking sequence used to make the repeat expansion region subject sequence (i.e 5000 for C9Orf72)
EXT_MIN = 10 # min number of GC bases beyond detected repeat region beyond which to consider extending 
MIN_EXT_FRAC = 0.95 # min GC fraction of the "flank" region needed to continue extension 
REP_UNIT_LEN = 6 # length of the repeat expansion unit

# parse reads 
record_dict = SeqIO.index(read_fasta_file, "fasta")

# parse the blast file, and store alignments in a defaultdict
@dataclass
class Aln:
	qstart: int
	qend: int
	sstart: int
	send: int
	score: float

read_alns = defaultdict(list)
with open(blast_outfmt6_file) as f:
	for line in f:
		tabs = line.rstrip('\n').split('\t')
		read_id = tabs[0]
		read_alns[read_id].append( Aln(int(tabs[6]),int(tabs[7]),int(tabs[8]),int(tabs[9]),float(tabs[11])) )

ALN_DIST_TO_FLANK_END = 100 # max dist from FLANK_LEN to consider an alignment as potentially bordering a repeat exp

def get_flank_and_rep_seqs(read_id,max_qcoords):
	lflank = record_dict[read_id].seq[max_qcoords[0]:max_qcoords[1]]
	rep_seq = record_dict[read_id].seq[max_qcoords[1]:max_qcoords[2]]
	rflank = record_dict[read_id].seq[max_qcoords[2]:max_qcoords[3]]
	if rep_seq.count('CCCCGG') > rep_seq.count('GGGGCC'):
		lflank = str( Seq(rflank).reverse_complement() )
		rep_seq = str( Seq(rep_seq).reverse_complement() )
		rflank = str( Seq(lflank).reverse_complement() )
	return (lflank,rep_seq,rflank)

def print_repeat_exp_seq(read_id,aln_list):
	if len(aln_list) < 2:
		return None
	## try to find a pair of alignments that corresponds to one with sstart near the FLANK_LEN, one with ssend near the FLANK_LEN, and both being in the same alignment ddirection with respect to the query
	sstarts = []
	sends = []
	for a in aln_list:
		if abs(a.sstart - FLANK_LEN) <= ALN_DIST_TO_FLANK_END:
			sstarts.append(a)
		elif abs(a.send - FLANK_LEN) <= ALN_DIST_TO_FLANK_END:
			sends.append(a)
	max_flank = 0
	max_qcoords = []
	for ss in sstarts:
		for sa in sends:
			if (ss.send - ss.sstart)*(sa.send - sa.sstart) > 0: # this checks if the two alignments are in the same direction (plus/minus) with respect to the query
				if abs(ss.send - ss.sstart) + abs(sa.send - sa.sstart) > max_flank:
					max_flank = abs(ss.send - ss.sstart) + abs(sa.send - sa.sstart)
					max_qcoords = sorted([ss.qstart,ss.qend,sa.qstart,sa.qend])
	if max_qcoords != []:
		lflank,rep_seq,rflank = get_flank_and_rep_seqs(read_id,max_qcoords)
		kd = Kmer_Dict(rep_seq)
		print(f"{read_id}\t{round(len(rep_seq)/REP_UNIT_LEN,2)}\t{lflank}\t{rep_seq}\t{rflank}\t{kd.get_condensed_seq()}")

for read_id in read_alns:
	print_repeat_exp_seq(read_id,read_alns[read_id])
