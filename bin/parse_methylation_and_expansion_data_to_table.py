#!/usr/bin/python3
import sys
from collections import defaultdict
import numpy as np
from sklearn.cluster import KMeans

expansion_file = sys.argv[1]
methylation_file = sys.argv[2] # format is modkit call-mods followed by extract 
lb = int(sys.argv[3]) # lower co-ordinate on chromosome to consider methylation information
ub = int(sys.argv[4]) # upper co-ordinate on chromosome to consider methylation information


def get_gc_perc(seq):
	if len(seq) == 0:
		return None
	seq = seq.upper()
	return round(100*(seq.count('G') + seq.count('C'))/len(seq),2)

## parse repeat expansions
read_exp_sizes = {}
with open(expansion_file) as f:
	for line in f:
		tabs = line.rstrip('\n').split('\t')
		if get_gc_perc(tabs[3]) < 75:
			continue
		read_exp_sizes[tabs[0]] = float(tabs[1])
	

## parse methylation info, keeping track of positions that are called in all reads
read_meth_posns = defaultdict(dict)
all_positions = set()
with open(methylation_file) as f:
	header = f.readline()
	for line in f:
		tabs = line.rstrip('\n').split('\t')
		if tabs[10] != 'm':
			continue
		pos = int(tabs[2])
		if (pos < lb) or (pos > ub):
			continue
		rid = tabs[0]
		meth = round(float(tabs[9]))
		if meth == 0:
			meth = -1
		read_meth_posns[rid][pos] = meth
		all_positions.add(pos)

read_count = len(read_meth_posns)
vecs_for_clustering = {}
for rid in sorted(read_exp_sizes, key = lambda x:read_exp_sizes[x]):
	outline = f"{rid}\t{read_exp_sizes[rid]}\t" 
	vector = []
	for pos in all_positions:
		output = 0
		if pos in read_meth_posns[rid]:
			output = read_meth_posns[rid][pos]
		outline += f"{output},"
		vector.append(output)
	#print(outline.rstrip(','))
	vecs_for_clustering[rid] = vector

def get_meth_frac(vec):
	num_sites = sum([1 for x in vec if x != 0])
	if num_sites == 0:
		return 0
	return sum([1 for x in vec if x == 1]) / num_sites

header = "read_id\texpansion_size\tmethylation_fraction\t" + "\t".join([str(pos) for pos in sorted(all_positions)])
print(header)
for rid in sorted(vecs_for_clustering, key = lambda x: get_meth_frac(vecs_for_clustering[x])):
	outline = f"{rid}\t{read_exp_sizes[rid]}\t{get_meth_frac(vecs_for_clustering[rid])}\t" + "\t".join([str(v) for v in vecs_for_clustering[rid]])
	print(outline)
