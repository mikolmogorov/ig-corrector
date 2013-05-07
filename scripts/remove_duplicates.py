#!/usr/bin/env python

import fasta_reader as fr
import sys

def remove_dups(seqs):
	counter = 0
	new_seqs = {}
	seq_count = {}
	for h in seqs:
		seq_count[seqs[h]] = seq_count.get(seqs[h], 0) + 1
	for s in seq_count:
		new_seqs["Seq{0}_{1}".format(counter, seq_count[s])] = s
		counter += 1
	return new_seqs

def main():
	if len(sys.argv) < 2:
		print "USAGE: remove_duplicates.py filename"
		return

	newSeqs = remove_dups(fr.get_seqs(sys.argv[1]))
	for h, seq in newSeqs.iteritems():
		print ">{0}\n{1}".format(h, seq)

if __name__ == "__main__":
	main()
