#!/usr/bin/env python

from random import sample, choice
import editdist
import subprocess
import sys
from math import log
import fasta_reader as fr
import msa


def k_means(sequences, k):
	centers = dict(enumerate(sample(sequences.values(), k)))
	old_centers = [[]]

	for i in xrange(0, 100):
		clusters = {i : [] for i in xrange(k)}
		sys.stderr.write(str(i) + " ")

		for header, seq in sequences.iteritems():
			min_dist = sys.maxint, 0
			for c_id, c_str in centers.iteritems():
				dist = get_distance(seq, c_str)
				if dist < min_dist[0]:
					min_dist = dist, c_id
			clusters[min_dist[1]].append(header)

		new_centers = {}
		for cl_id, cl in clusters.iteritems():
			if len(cl) > 0:
				new_centers[cl_id] = msa.get_consensus(cl, sequences)
			else:
				new_centers[cl_id] = choice(sequences.values())

		shift = 0
		for i in xrange(k):
			if centers[i] != new_centers[i]:
				shift += 1
		if shift == 0:
			break

		#print "shifted centers: ", shift
		#if new_centers in old_centers:
			#print "Stuck in a loop, trying new guess"
		#	return k_means(sequences, k)

		centers = new_centers
		old_centers.append(centers)
	sys.stderr.write("\n")
	return clusters


def get_distance(seq1, seq2):
	key = seq1 + seq2 if seq1 < seq2 else seq2 + seq1
	if key not in get_distance.dist_cache:
		get_distance.dist_cache[key] = editdist.distance(seq1, seq2)
	return get_distance.dist_cache[key]
get_distance.dist_cache = {}


def main():
	assert len(sys.argv) >= 3
	seqs = fr.get_seqs(sys.argv[1])

	k = int(sys.argv[2])
	clusters = k_means(seqs, k)
	for cl in clusters.values():
		for h in cl:
			print ">{0}\n{1}".format(h, seqs[h])
		print ""


if __name__ == "__main__":
	main()
