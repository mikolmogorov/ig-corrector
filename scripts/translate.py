#!/usr/bin/env python

from Bio.Seq import Seq
import fasta_reader as fr
import sys

TRHLD = 1
for h, seq in fr.read_fasta(open(sys.argv[1], "r")).iteritems():
	mult = int(h.split("_")[2]) if len(h.split("_")) > 2 else 1
	if mult >= TRHLD:
		print ">{0}\n{1}".format(h, str(Seq(seq[int(sys.argv[2]):]).translate()))
