import sys
from collections import defaultdict
from dataclasses import dataclass

class Kmer_Dict:
	"""Methods to create a kmer dictionary of the input string and find consecutively repeated kmers"""

	def __init__(self,seq,KMER_SIZE=6,MIN_CONDENSE_REPS=3):
		self.seq = seq
		self.KMER_SIZE=KMER_SIZE
		self.MIN_CONDENSE_REPS=MIN_CONDENSE_REPS
		self.kdict = self._build_dict() # dict of kmer sequence to list of starting indices
		self._prune_kdict_by_permutation_groups()

	def _build_dict(self):
		kdict = defaultdict(list)
		for i in range(len(self.seq)):
			for j in range(i+self.KMER_SIZE,i+self.KMER_SIZE+1):
				kdict[ self.seq[i:j] ].append(i) 
		return kdict

	def _get_permutation_set_from_kmer(self,kmer):
		pgroup = set()
		kmer2 = kmer + kmer
		for i in range(len(kmer)):
			pgroup.add(kmer2[i:i + len(kmer)])
		return pgroup

	def _prune_kdict_by_permutation_groups(self):
		pruned_kdict = {}
		kmer_pgroups = []
		for kmer in sorted(self.kdict,key = lambda x: len(self.kdict[x])):
			kmer_pgroups.append(self._get_permutation_set_from_kmer(kmer))
		for pg in kmer_pgroups:
			max_occ = 0
			max_kmer = None
			for kmer in pg:
				if len(self.kdict[kmer]) > max_occ:
					max_occ = len(self.kdict[kmer])
					max_kmer = kmer
			pruned_kdict[max_kmer] = self.kdict[max_kmer]
		self.kdict = pruned_kdict

	def _get_repeated_kmer_position_dict(self):
		"""Return a dict of position: [number of repeats] for all kmers that repeat at least self.MIN_CONDENSE_REPS times"""
		pos_dict = {}
		for kmer in self.kdict:
			in_rep = False
			num_reps = 1
			rep_start_pos = 0
			for i in range(len(self.kdict[kmer])-1):
				if self.seq[self.kdict[kmer][i] : self.kdict[kmer][i+1]] == kmer:
					num_reps += 1
					if not in_rep:
						in_rep = True
						rep_start_pos = self.kdict[kmer][i]
				elif in_rep:
					if num_reps >= self.MIN_CONDENSE_REPS:
						pos_dict[rep_start_pos] = num_reps
					in_rep = False
					num_reps = 1
					rep_start_pos = 0
			## close the final rep
			if in_rep and (num_reps >= self.MIN_CONDENSE_REPS):
				pos_dict[rep_start_pos] = num_reps
		return pos_dict

	def get_condensed_seq(self):
		rep_pos_dict = self._get_repeated_kmer_position_dict()
		if len(rep_pos_dict) == 0:
			return self.seq
		rep_pos = [pos for pos in sorted(rep_pos_dict)]
		num_rep_list = [rep_pos_dict[pos] for pos in sorted(rep_pos_dict)]
		condensed_seq	= self.seq[:rep_pos[0]]
		for i in range(len(rep_pos)-1):
			num_reps = num_rep_list[i]
			condensed_seq += "(" + self.seq[ rep_pos[i] : rep_pos[i] + self.KMER_SIZE ] + f"[{num_reps}])"
			condensed_seq += self.seq[ rep_pos[i] + self.KMER_SIZE*num_reps : rep_pos[i+1] ]
		condensed_seq += "(" + self.seq[ rep_pos[-1] : rep_pos[-1] + self.KMER_SIZE ] + f"[{num_rep_list[-1]}])"
		condensed_seq += self.seq[rep_pos[-1] + self.KMER_SIZE*num_rep_list[-1] : ]
		return condensed_seq
									

if __name__ == "__main__":
	seq = sys.argv[1]
	kd = Kmer_Dict(seq)
	#for kmer in sorted(kd.kdict,key = lambda x: len(kd.kdict[x])):
	#	print(kmer,kd.kdict[kmer])
		#print(kd._get_permutation_set_from_kmer(kmer))
	print("\n")
	print(kd.get_condensed_seq())

