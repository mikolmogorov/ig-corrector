#!/usr/bin/env python

import fasta_reader as fr
import msa
import sys
import editdist
from Bio import pairwise2
from itertools import product, combinations
from collections import defaultdict

def merge_cdrs(cdr_map, weight, full_seqs, threshold):
	cons_cache = {}
	sorted_r = sorted(weight, key = weight.get, reverse = True)
	for cdr1, cdr2 in combinations(sorted_r, 2):
		if cdr1 not in cdr_map or cdr2 not in cdr_map:
			continue

		align = pairwise2.align.globalms(cdr1, cdr2, 0, -1, -2, -2)[0]
		n_miss = 0
		for i in xrange(len(align[0])):
			if align[0][i] != align[1][i] and align[0][i] != "-" and align[1][i] != "-":
				n_miss += 1

		#sys.stderr.write(".")
		if n_miss > 0:
			if cdr1 not in cons_cache:
				cons_cache[cdr1] = msa.get_consensus(cdr_map[cdr1], full_seqs)
			if cdr2 not in cons_cache:
				cons_cache[cdr2] = msa.get_consensus(cdr_map[cdr2], full_seqs)
			dist = editdist.distance(cons_cache[cdr1], cons_cache[cdr2])
			#sys.stderr.write(" " + str(dist) + " ")
			if dist > 4:
				continue
		
		if float(weight[cdr1]) / weight[cdr2] >= threshold:
			true_cdr, false_cdr = cdr1, cdr2
		elif len(cdr1) % 3 == 0 and len(cdr2) % 3 != 0:
			true_cdr, false_cdr = cdr1, cdr2
		elif len(cdr2) % 3 == 0 and len(cdr1) % 3 != 0:
			true_cdr, false_cdr = cdr2, cdr1
		else:
			sys.stderr.write("Ambigious cdrs! " + str(align) + "\n")
			true_cdr, false_cdr = cdr1, cdr2
		#print false_cdr, "to", true_cdr, weight[false_cdr], weight[true_cdr], len(cdr_map)

		cdr_map[true_cdr].extend(cdr_map[false_cdr])
		del cdr_map[false_cdr]
		weight[true_cdr] += weight[false_cdr]
		del weight[false_cdr]
		cons_cache[true_cdr] = msa.get_consensus(cdr_map[true_cdr], full_seqs)


def correct_cluster(clust_seqs, full_seqs, threshold):
	#cdr_map = {}
	cdr_map = defaultdict(list)
	weight = {}
	for seq_head, seq_cdr in clust_seqs.iteritems():
		#if not seq_cdr in cdr_map:
		#	cdr_map[seq_cdr] = []
		cdr_map[seq_cdr].append(seq_head)
		qty = int(seq_head.split("_")[1])
		weight[seq_cdr] = weight.get(seq_cdr, 0) + qty
	
	if len(cdr_map) == 1:
		return {"Cluster" : clust_seqs}
	
	merge_cdrs(cdr_map, weight, full_seqs, threshold)

	clusters = {}
	counter = 0
	for cdr in cdr_map:
		clust_seqs = {}
		for head in cdr_map[cdr]:
			clust_seqs[head] = cdr
		clusters["Cluster{0}".format(counter)] = clust_seqs
		counter += 1
	return clusters
	

def out_cluster(name, clust_seqs, full_seqs, out_stream):
	out_stream.write("=" + name + "\n")
	for hdr, seq in clust_seqs.iteritems():
		out_stream.write(">{0} {1}\n{2}\n".format(hdr, seq, full_seqs[hdr]))


def correct_cdr(cdr_stream, seqs_stream, threshold, out_stream):
	full_seqs = fr.read_fasta(seqs_stream)
	clusters = fr.read_cluster(cdr_stream)
	cl_id = 0
	for c_name, c_seqs in clusters.iteritems():
		for newc_name, newc_seqs in correct_cluster(c_seqs, full_seqs, threshold).iteritems():
			out_cluster("Cluster_{0}".format(cl_id), newc_seqs, full_seqs, out_stream)
			cl_id += 1


def main():
	THRESHOLD = 2
	if len(sys.argv) < 3:
		print "USAGE: correct_cdr.py cdr_cluster_file reads_file"
		return
	correct_cdr(open(sys.argv[1], "r"), open(sys.argv[2], "r"), THRESHOLD, sys.stdout)
	

if __name__ == "__main__":
	main()
