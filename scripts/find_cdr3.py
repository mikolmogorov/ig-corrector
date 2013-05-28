#!/usr/bin/env python

import dna_vs_prot as aln
import fasta_reader as fr
import sys
from Bio.Seq import Seq


CDR3_START = "YYC"
CDR3_END = "WG[QKRS]"
#CDR3_END = "GTK"


def main():
	graph_start = aln.build_graph(CDR3_START)
	graph_end = aln.build_graph(CDR3_END)

	counter = 0

	for h, seq in fr.get_seqs(sys.argv[1]).iteritems():
		counter += 1
		cdr_start = aln.loc_align(seq, graph_start, 16)
		cdr_end = aln.loc_align(seq, graph_end, 16)
		candidates = []
		for pos_start in cdr_start:
			for pos_end in cdr_end:
				dist = pos_end[0] - pos_start[0]
				if 20 <= dist and dist <= 80 and pos_start[0] > len(seq) / 2:
					candidates.append((pos_start[0], pos_end[1], pos_start[2] + pos_end[2]))

		#print seq
		if len(candidates) == 0:
			#print "cdr3 not found"
			continue
		
		maxScore = 0
		cand = None
		for c in candidates:
			if c[2] > maxScore:
				cand = c
				maxScore = c[2]
		#for c in candidates:
		cdr = seq[cand[0] : cand[1] + 1]
		sys.stdout.write(">{0}\n{1}\n".format(h, cdr))
		sys.stderr.write(str(counter) + "\n")
		#print counter, cand, cdr, str(Seq(cdr).translate())


if __name__ == "__main__":
	main()
