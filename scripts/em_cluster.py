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
	divergence = max(pow(eps, edit), 0.00000000001)
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

	sigma = 1.0 / 60

	#while True:
	for _ in xrange(100):
		#print "step"
		#M step
		U = [0] * n_seqs
		for j in xrange(n_seqs):
			s_min = sys.maxint, 0
			for k in xrange(n_seqs):
				s_dist = 0
				for i in xrange(n_seqs):
					s_dist += z[i][j] * distances[i][k]#get_distance(sequences[i], sequences[k])
				if s_dist < s_min[0]:
					s_min = s_dist, k
			#U[j] = sequences[s_min[1]]
			U[j] = s_min[1]

		tau = [0.0] * n_seqs
		for i in xrange(n_seqs):
			for j in xrange(n_seqs):
				tau[i] += z[j][i]
			tau[i] /= n_seqs
		#print tau

		#E step
		for i in xrange(n_seqs):
			for j in xrange(n_seqs):
				#up = tau[j] * exp(-get_distance(sequences[i], U[j]) / sigma)
				up = tau[j] * exp(-distances[i][U[j]] / sigma)
				norm = 0.0
				for k in xrange(n_seqs):
					#norm += tau[k] * exp(-get_distance(sequences[i], U[k]) / sigma)
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
		clusters[p.n].append(seq_enum[n.n])
	
	for c in clusters:
		fr.write_fasta(dict(clusters[c]), sys.stdout)
		print ""
	
	#for y in xrange(n_seqs):
	#	print "%f" % z[x][y],
	#print ""

def main():
	seqs = fr.get_seqs(sys.argv[1])
	em_cluster(seqs)


if __name__ == "__main__":
	main()
