#!/usr/bin/env python2

import sys
from Bio.Seq import Seq
from Bio import pairwise2

def fastq_source(filename):
	fin = open(filename, "r")
	while True:
		header = fin.readline().strip()
		if not header:
			break
		seq = fin.readline().strip()
		fin.readline()
		qual = fin.readline().strip()
		yield header, seq, qual


mids = {"VH" : ("ATAGATAGAC", "ATATAGTCGC"), 
		"VL" : ("ATCTACTGAC", "CACGTAGATC"),
		"VK" : ("CACGTGTCGC", "CATACTCTAC")}

primers = {	"VH" : ("GCCTACGGCAGCCGCTGGATTGTTATTAC", "CACAGACGGGCCTTTTGTAGAC"),
			"VL" : ("GCTGCTGCTGGTCTGCTGCTCCTCGCTG", "GGCGGGAAAATAAAAACAGACGG"),
			"VK" : ("GCTGCTGCTGGTCTGCTGCTCCTCGCTG", "GGCGGGAAAATAAAAACAGACGG")}


def with_sequence(m, seq):
	rev_seq = str(Seq(seq).reverse_complement())

	for i in [0, 1]:
		start_seq = mids[m][i] + primers[m][i]
		start_aln = pairwise2.align.localms(seq[0:50], start_seq, 1, -1, -1, -1)[0]
		if start_aln[2] < len(start_seq) - 3:
			#print "fail start"
			continue
		start = len(start_aln[1].rstrip("-"))

		end_seq = mids[m][abs(1 - i)] + primers[m][abs(1 - i)]
		end_aln = pairwise2.align.localms(rev_seq[0:200], end_seq, 1, -1, -1, -1)[0]
		if end_aln[2] < len(end_seq) - 9:
			#print "fail end"
			continue
		end = len(end_aln[1].rstrip("-"))
		
		res = seq[start : -end]
		if i == 0:
			return res
		else:
			return str(Seq(res).reverse_complement())

	return None



def main():
	processed = 0
	files = {}
	for mid in mids:
		files[mid] = open(mid + ".fasta", "w")

	files["None"] = open("no-mid.fastq", "w")

	for header, seq, qual in fastq_source(sys.argv[1]):
		print processed
		processed += 1

		count = 0

		for m in mids:
			str_seq = seq.strip("acgtn")
			res = with_sequence(m, str_seq)
			if res:
				files[m].write("{0}\n{1}\n".format(">" + header, res))
				count += 1
		
		assert count <= 1
		if count == 0:
			files["None"].write("{0}\n{1}\n+\n{2}\n".format(header, seq, qual))

if __name__ == "__main__":
	main()
