#!/usr/bin/env python2

import editdist
import sys
import fasta_reader as fr
from math import log, exp

##disjoint set
class Node:
	def __init__(self, nn):
		self.n = nn
		self.parent = self
		self.rank = 0

def make_set(x):
      x.parent = x
      x.rank   = 0
 
def union(x, y):
	x_root = find(x)
	y_root = find(y)
	if x_root == y_root:
		return
 
	if x_root.rank < y_root.rank:
		x_root.parent = y_root
	elif x_root.rank > y_root.rank:
		y_root.parent = x_root
	else:
		y_root.parent = x_root
		x_root.rank = x_root.rank + 1

def find(x):
	if x.parent != x:
		x.parent = find(x.parent)
	return x.parent
##############


def get_distance(seq1, seq2):
	eps = 0.01 #magic
	edit = editdist.distance(seq1, seq2)
	avg_len = float(len(seq1) + len(seq2)) / 2
	divergence = max(pow(eps, edit), sys.float_info.epsilon)
	dist = log(pow(1 - eps, avg_len - edit)) + log(divergence)
	return -dist / avg_len

def init_z(sequences):
	length = len(sequences)
	matrix = [[0 if x != y else 1 for x in xrange(length)] for y in xrange(length)] 
	return matrix

def precalc_dist(sequences):
	length = len(sequences)
	matrix = [[get_distance(sequences[x], sequences[y]) for x in xrange(length)] for y in xrange(length)]
	return matrix

def em_cluster(fasta_dict):
	seq_enum = {i : s for i, s in enumerate(fasta_dict.iteritems())}
	sequences = [seq_enum[i][1] for i in seq_enum]

	z = init_z(sequences)
	n_seqs = len(sequences)
	distances = precalc_dist(sequences)

	sigma = 1.0 / 80

	#while True:
	for _ in xrange(20):
		#print "step"
		#M step
		U = [0] * n_seqs
		for j in xrange(n_seqs):
			s_min = sys.maxint, 0
			for k in xrange(n_seqs):
				s_dist = 0
				for i in xrange(n_seqs):
					s_dist += z[i][j] * distances[i][k]
				if s_dist < s_min[0]:
					s_min = s_dist, k
			U[j] = s_min[1]

		tau = [0.0] * n_seqs
		for i in xrange(n_seqs):
			for j in xrange(n_seqs):
				tau[i] += z[j][i]
			tau[i] /= n_seqs

		#E step
		for i in xrange(n_seqs):
			for j in xrange(n_seqs):
				up = tau[j] * exp(-distances[i][U[j]] / sigma)
				norm = 0.0
				for k in xrange(n_seqs):
					norm += tau[k] * exp(-distances[i][U[k]] / sigma)
				z[i][j] = up / norm

	origins = {}
	nodes = [Node(x) for x in xrange(n_seqs)]
	for x in xrange(n_seqs):
		origins[x] = max( (v,i) for i, v in enumerate(z[x]) )[1]
		union(nodes[x], nodes[origins[x]])

	clusters = {}
	for n in nodes:
		p = find(n)
		if p.n not in clusters:
			clusters[p.n] = []
		clusters[p.n].append(seq_enum[n.n][0])
	
	return clusters
	

def main():
	seqs = fr.get_seqs(sys.argv[1])
	clusters = em_cluster(seqs)

	for cl in clusters.values():
		for h in cl:
			print ">{0}\n{1}".format(h, seqs[h])
		print ""

	#print len(clusters)

if __name__ == "__main__":
	main()
