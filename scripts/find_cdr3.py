#!/usr/bin/env python

import dna_vs_prot as aln
import fasta_reader as fr
import sys
from Bio.Seq import Seq

CDR3_START = "YYC"
CDR3_END = "W[GS][QKRS]"

graph_start = aln.build_graph(CDR3_START)
graph_end = aln.build_graph(CDR3_END)

for h, seq in fr.fasta_source(sys.argv[1]):
	cdr_start = aln.loc_align(seq, graph_start, 16)
	cdr_end = aln.loc_align(seq, graph_end, 16)
	candidates = []
	for pos_start in cdr_start:
		for pos_end in cdr_end:
			dist = pos_end[0] - pos_start[0]
			if 20 <= dist and dist <= 80 and pos_start[0] > len(seq) / 2:
				candidates.append((pos_start[0], pos_end[1], pos_start[2] + pos_end[2]))

	#print seq
	#assert len(candidates) <= 1
	if len(candidates) == 0:
		print "cdr3 not found"
		continue
	
	maxScore = 0
	cand = None
	for c in candidates:
		if c[2] > maxScore:
			cand = c
			maxScore = c[2]
	#for c in candidates:
	cdr = seq[cand[0] : cand[1] + 1]
	print cand, cdr, str(Seq(cdr).translate())
