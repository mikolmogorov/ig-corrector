#!/usr/bin/env python

import sys
import fasta_reader as fr

#s = set()
#for h, _ in fr.fasta_source(sys.argv[1]):
#	assert h not in s
#	s.add(h)

#print len(s)

for line in open(sys.argv[1], "r"):
	if line.startswith("Cluster") and int(line.split("_")[1]) > 1:
		print line.strip("\n")
