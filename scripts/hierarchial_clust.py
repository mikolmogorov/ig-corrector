#!/usr/bin/env python

import fasta_reader as fr
import msa
import sys
import editdist
from itertools import product, combinations

class Cluster:
	def __init__(self):
		self.seqs = []

def cluster_dist(clust1, clust2, seqs):
	dist_sum = 0.0
	for h1, h2 in product(clust1.seqs, clust2.seqs):
		#print seq1, seq2
		dist_sum += seq_disstance(h1, h2, seqs)
		#print dist_sum
	return dist_sum / (len(clust1.seqs) * len(clust2.seqs))


def seq_disstance(seq1, seq2, seqs):
	key = dict_key(seq1, seq2)
	if key not in seq_disstance.cache:
		seq_disstance.cache[key] = editdist.distance(seqs[seq1], seqs[seq2])
	return seq_disstance.cache[key]
seq_disstance.cache = {}


def cons_score(heads, seqs):
	align = msa.align_muscle(heads, seqs)
	score = 0
	seq_len = len(align[align.keys()[0]])
	for i in xrange(seq_len):
		sym = align[align.keys()[0]][i]
		fail = False
		for h in align:
			if align[h][i] != sym:
				fail = True
				break
		if not fail:
			score += 1
	return float(score) / seq_len * len(heads)
			

def calc_score(clusters, seqs):
	score = 0
	for cl in clusters:
		score += len(clusters[cl].seqs) * cons_score(clusters[cl].seqs, seqs) 
	return score


def dict_key(c1, c2): 
	return tuple(sorted([c1, c2]))


def step(distances, clusters, seqs, cutoff):
	cl1, cl2 = min(distances, key = distances.get)
	key = dict_key(cl1, cl2)

	if distances[key] > cutoff:
		return False

	sys.stderr.write("merging {0} {1} {2}\n".format(cl1, cl2, distances[key]))
	assert cl1 in clusters and cl2 in clusters

	clusters[cl1].seqs += clusters[cl2].seqs
	for cl in clusters:
		#print dict_key(cl, cl2)
		if cl != cl2:
			del distances[dict_key(cl, cl2)]
		if cl != cl1 and cl != cl2:
			distances[dict_key(cl, cl1)] = cluster_dist(clusters[cl], clusters[cl1], seqs)
		assert dict_key(cl, cl2) not in distances
	del clusters[cl2]

	return True


def cluster(seqs, cutoff):
	clusters = {}
	counter = 0
	for s in seqs:
		clusters[counter] = Cluster()
		clusters[counter].seqs.append(s)
		counter += 1
	
	distances = {}
	for cl1, cl2 in combinations(clusters, 2):
		d = cluster_dist(clusters[cl1], clusters[cl2], seqs)
		distances[dict_key(cl1, cl2)] = d

	while True:
		res = step(distances, clusters, seqs, cutoff)
		if not res or len(clusters) == 1:
			break
		#sys.stderr.write(str(calc_score(clusters, seqs)) + "\n")
	
	return {n : cl.seqs for n, cl in clusters.iteritems()}
	

def main():
	if len(sys.argv) < 3:
		sys.stderr.write("USAGE: hierarchial_clust.py fasta_file cutoff\n")
		return -1
		
	seqs = fr.get_seqs(sys.argv[1])
	clusters = cluster(seqs, float(sys.argv[2]))

	for cl in clusters.values():
		for h in cl:
			print ">{0}\n{1}".format(h, seqs[h])
		print ""


if __name__ == "__main__":
	main()
