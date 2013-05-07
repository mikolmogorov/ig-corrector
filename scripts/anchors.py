#!/usr/bin/env python2

import fasta_reader as fr
import sys
from itertools import combinations

KMER_LEN = 21
THRESHOLD = 0.5

#disjoint sets
class Node:
	def __init__(self, data):
		self.parent = self
		self.rank = 0
		self.data = data

def make_set(x):
	return Node(x)

def union(x, y):
     xRoot = find(x)
     yRoot = find(y)
     if xRoot.rank > yRoot.rank:
         yRoot.parent = xRoot
     elif xRoot.rank < yRoot.rank:
         xRoot.parent = yRoot
     elif xRoot != yRoot: # Unless x and y are already in same set, merge them
         yRoot.parent = xRoot
         xRoot.rank = xRoot.rank + 1

def find(x):
     if x.parent == x:
        return x
     else:
        x.parent = find(x.parent)
        return x.parent
##########

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

def extract_kmers(seqs, kmer_len):
	kmer_hash = {}
	seqs_kmers = {}
	counter = 0
	for header in seqs:
		seqs_kmers[header] = []
		for kmer in fr.iter_kmers(seqs[header], kmer_len):
			kmer_id = None
			if kmer in kmer_hash:
				kmer_id = kmer_hash[kmer]
			else:
				kmer_hash[kmer] = counter
				kmer_id = counter
				counter += 1
			seqs_kmers[header].append(kmer_id)
	#for h in seqs_kmers:
	#	seqs_kmers[h] = sorted(seqs_kmers[h])

	return seqs_kmers

def cluster_seqs(seqs, kmer_len, threshold):
	kmers = extract_kmers(seqs, kmer_len)
	#avg_sim = 0
	n_seqs = 0
	#f = open("strange.fasta", "w")
	nodes = { c : make_set(c) for c in seqs }
	old_precent = 0

	for seq1, seq2 in combinations(seqs, 2):
		first_kmers = set()
		for kmer in kmers[seq1]:
			first_kmers.add(kmer)
		counter = 0
		for kmer in kmers[seq2]:
			if kmer in first_kmers:
				counter += 1
		
		n_kmers = (len(kmers[seq1]) + len(kmers[seq2])) / 2
		sim_rate = float(counter) / n_kmers
		#avg_sim += sim_rate
		n_seqs += 1

		if sim_rate > threshold:
			union(nodes[seq1], nodes[seq2])

		percent = 100 * n_seqs / (len(seqs) * len(seqs) / 2)
		if percent > old_precent:
			print percent
			old_precent = percent

		#print sim_rate
		#if counter == 0:
		#	f.write(seqs[seq1] + "\n" + seqs[seq2] + "\n\n")	
		#print n_seqs, " %6.2f" % (float(n_seqs) * 100 / (len(seqs) * len(seqs) / 2))

	clusters = {}
	for seq in nodes:
		cur_set = find(nodes[seq]).data
		if not cur_set in clusters:
			clusters[cur_set] = []
		clusters[cur_set].append(seq)

	return clusters
	#print avg_sim / n_seqs, len(seqs)
	#f.close()

def test():
	seqs = fr.get_seqs("alone.fasta")
	headers = []
	cur_align = None
	for line in open("aligns.txt", "r"):
		line = line.strip("\n")
		if len(line) == 0:
			if len(headers) > 0:
				headers.append(cur_align)
				clust_seqs = {h : seqs[h] for h in headers}

				cluster_seqs(clust_seqs, KMER_LEN)
			headers = []
		elif line[0] == ">":
			cur_align = line[1:]
		else:
			headers.append(line.split(" ")[0])	

def test2():
	seqs = fr.get_seqs("MID86-F.fasta")
	trimmed = remove_dups(seqs)
	clust = cluster_seqs(trimmed, KMER_LEN, THRESHOLD)
	out_file = open("clusters.txt", "w")
	for i, cl in enumerate(clust):
		out_file.write("Cluster{0}_{1}\n".format(i, len(clust[cl])))
		for seq in clust[cl]:
			out_file.write(">{0}\n{1}\n".format(seq, trimmed[seq]))
		out_file.write("\n")

if __name__ == "__main__":
	#s = make_set("aa")
	#ss = make_set("bb")
	#union(s, ss)

	#print find(ss).data
	test2()
