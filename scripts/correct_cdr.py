#!/usr/bin/env python

import fasta_reader as fr
import sys
from Bio import pairwise2
from itertools import product, combinations

class Sequence:
	def __init__(self, hdr, cdr):
		self.header = hdr
		self.cdr = cdr

class Cluster:
	def __init__(self):
		self.seqs = []

def parse_cdr(filename):
	clusters = []
	cl = None
	header = ""
	for line in open(filename, "r"):
		line = line.strip()
		if line.startswith("Cluster"):
			if cl:
				clusters.append(cl)
			cl = Cluster()
		elif line.startswith(">"):
			header = line[1:]
		elif len(line) > 0:
			cl.seqs.append(Sequence(header, line))
	return clusters

def correct_cdr(cluster, seqs):
	THRESHOLD = 2
	cdr_map = {}
	weight = {}
	for s in cluster.seqs:
		if not s.cdr in cdr_map:
			cdr_map[s.cdr] = []
		cdr_map[s.cdr].append(s.header)
		qty = int(s.header.split("_")[1])
		weight[s.cdr] = weight.get(s.cdr, 0) + qty
	
	#print cluster.seqs
	if len(cdr_map) == 1:
		return cluster
	print len(cdr_map)

	#correct indels
	done = False
	while not done:
		#sorted_f = sorted(weight, key = weight.get)
		sorted_r = sorted(weight, key = weight.get, reverse = True)
		done = True

		#for cdr1, cdr2 in zip(sorted_r, sorted_f):
		#for cdr1, cdr2 in product(sorted_r, sorted_f):
		for cdr1, cdr2 in combinations(sorted_r, 2):
			if cdr1 == cdr2:
				continue

			align = pairwise2.align.globalms(cdr1, cdr2, 0, -1, -2, -2)[0]
			#print align
			n_miss = 0
			for i in xrange(len(align[0])):
				if align[0][i] != align[1][i] and align[0][i] != "-" and align[1][i] != "-":
					n_miss += 1
			if n_miss != 0:
				continue
			
			done = False
			
			if weight[cdr1] / weight[cdr2] >= THRESHOLD:
				true_cdr, false_cdr = cdr1, cdr2
			elif len(cdr1) % 3 == 0 and len(cdr2) % 3 != 0:
				true_cdr, false_cdr = cdr1, cdr2
			elif len(cdr2) % 3 == 0 and len(cdr1) % 3 != 0:
				true_cdr, false_cdr = cdr2, cdr1
			else:
				print "Ambigious cdrs!"
				print align
				true_cdr, false_cdr = cdr1, cdr2
			print false_cdr, "to", true_cdr, weight[false_cdr], weight[true_cdr], len(cdr_map)

			cdr_map[true_cdr] += cdr_map[false_cdr]
			del cdr_map[false_cdr]
			weight[true_cdr] += weight[false_cdr]
			del weight[false_cdr]
			break

def out_cluster():
	seqs = fr.get_seqs(sys.argv[2])
	for line in open(sys.argv[1], "r"):
		line = line.strip()
		if line.startswith("Cluster"):
			sys.stdout.write("\n" + line + "\n")
		elif line.startswith(">"):
			header = line[1:]
		elif len(line) > 0:
			sys.stdout.write(">{0} {1}\n{2}\n".format(header, line, seqs[header]))
		#if line.startswith("Cluster") and int(line.split("_")[1]) > 1:
			#print line.strip("\n")	

def main():
	if len(sys.argv) < 3:
		print "USAGE: postproc.py cdr_cluster_file reads_file"
		return

	seqs = fr.get_seqs(sys.argv[2])
	clusters = parse_cdr(sys.argv[1])
	counter = 0
	for c in clusters:
		print counter
		counter += 1
		correct_cdr(c, seqs)
	
	
if __name__ == "__main__":
	main()
