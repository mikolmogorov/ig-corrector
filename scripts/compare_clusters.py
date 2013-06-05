#!/usr/bin/env python

import k_means_cluster as km
import em_cluster as em
import hierarchial_clust as hc
import sys
import editdist
from itertools import combinations


def dict_key(c1, c2): 
	return tuple(sorted([c1, c2]))


def seq_disstance(seq1, seq2, seqs):
	key = dict_key(seq1, seq2)
	if key not in seq_disstance.cache:
		seq_disstance.cache[key] = editdist.distance(seqs[seq1], seqs[seq2])
	return seq_disstance.cache[key]
seq_disstance.cache = {}


def parse_cluster(stream):
	clusters = []
	header = ""
	fasta_dict = None
	for line in stream:
		line = line.strip()
		if line.startswith("Cluster"):
			if fasta_dict:
				clusters.append(fasta_dict)
			fasta_dict = {}
		elif line.startswith(">"):
			header = line[1:].split(" ")[0]
		elif len(line) > 0:
			fasta_dict[header] = line
	return clusters


def cluster_score(cluster, seqs):
	score = 0.0
	for cl in cluster:
		avg_dist = 0.0
		for seq1, seq2 in combinations(cluster[cl], 2):
			avg_dist += seq_disstance(seq1, seq2, seqs)
		score += avg_dist / len(cluster[cl])
	return score


def main():
	preclusters = parse_cluster(open(sys.argv[1], "r"))
	for fasta_dict in preclusters:
		if len(fasta_dict) > 100 or len(fasta_dict) <= 3:
			continue
		print "size", len(fasta_dict)
		em_cluster = em.em_cluster(fasta_dict)
		hier_cluster = hc.cluster(fasta_dict, 4.0)
		kmeans_cluster = km.k_means(fasta_dict, len(hier_cluster))
		print "EM: {0}\tHIER: {1}".format(len(em_cluster), len(hier_cluster))
		print "EM: {0}\tHIER: {1}\tKM: {2}".format(cluster_score(em_cluster, fasta_dict),
												cluster_score(hier_cluster, fasta_dict),
												cluster_score(kmeans_cluster, fasta_dict))
		print ""


if __name__ == "__main__":
	main()
