#!/usr/bin/env python

import sys
import subprocess
import fasta_reader as fr
import msa

class Sequence:
	def __init__(self, hdr, seq):
		self.header = hdr
		self.seq = seq

class Cluster:
	def __init__(self):
		self.seqs = []

def parse_cluster(stream):
	clusters = []
	cl = None
	header = ""
	for line in stream:
		line = line.strip()
		if line.startswith("Cluster"):
			if cl:
				clusters.append(cl)
			cl = Cluster()
		elif line.startswith(">"):
			header = line[1:].split(" ")[0]
		elif len(line) > 0:
			cl.seqs.append(Sequence(header, line))
	return clusters


def split_cluster(cluster, seqs):
	#return [cluster]
	fasta_dict = {s.header : s.seq for s in cluster.seqs}
	cmdline = ["../src/fasta_clusta", "-k", "21", "-m", "4"]
	child = subprocess.Popen(cmdline, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
	fr.write_fasta(fasta_dict, child.stdin)
	#fr.write_fasta(fasta_dict, open("dump.fasta", "w"))
	child.stdin.close()
	#for line in child.stderr:
	#	sys.stderr.write(line)
	out_clusters = parse_cluster(child.stdout)
	return out_clusters


def main():
	seqs = fr.read_fasta(open(sys.argv[2], "r"))
	out_seqs = {}
	count = 0
	corrected_seqs = {}
	init_clusters = parse_cluster(open(sys.argv[1], "r"))
	for cl in init_clusters:
		sys.stderr.write(str(len(cl.seqs)) + " ")
		clusters = split_cluster(cl, seqs)
		for c in clusters:
			heads = [s.header for s in c.seqs]
			cons = msa.get_consensus(heads, seqs)
			size = 0
			for s in c.seqs:
				size += int(s.header.split("_")[1]) 
			out_seqs["Seq_{0}_{1}".format(count, size)] = cons
			count += 1
			sys.stderr.write(str(count) + "\n")
	
	fr.write_fasta(out_seqs, sys.stdout)


if __name__ == "__main__":
	#seqs = fr.read_fasta(open("dump-muscle.fasta", "r"))
	#heads = seqs.keys()
	#cons = msa.get_consensus(heads, seqs)
	main()
