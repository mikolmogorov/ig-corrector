#!/usr/bin/env python

from Bio.Seq import Seq
import fasta_reader as fr
import sys

for h, seq in fr.get_seqs(sys.argv[1]).iteritems():
	if int(h.split("_")[2]) > 0:
		print ">{0}\n{1}".format(h, str(Seq(seq[2:]).translate()))
