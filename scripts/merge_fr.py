#!/usr/bin/env python

import fasta_reader as fr
from Bio.Seq import Seq
import sys

def main():
	if len(sys.argv) < 3:
		print "USAGE: remove_duplicates.py forward_reads reverse_reads"
		return

	forward = fr.get_seqs(sys.argv[1])
	reverse = fr.get_seqs(sys.argv[2])
	for h, seq in forward.iteritems():
		print ">{0}\n{1}".format(h, seq)
	for h, seq in reverse.iteritems():
		print ">{0}\n{1}".format(h, str(Seq(seq).reverse_complement()))

if __name__ == "__main__":
	main()
