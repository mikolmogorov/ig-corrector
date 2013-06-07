#!/usr/bin/env python

import sys
import subprocess
import fasta_reader as fr
import msa


GRAPH_PATH = "graph_clust"
HIERARCH_PATH = "hierarchial_clust"


def run_graph(cluster_seqs, threshlod):
	cmdline = [GRAPH_PATH, "-k", "21", "-m", str(threshlod)]
	child = subprocess.Popen(cmdline, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
	fr.write_fasta(cluster_seqs, child.stdin)
	child.stdin.close()
	preclusters = fr.read_cluster(child.stdout)
	return preclusters


def run_hierarchial(cluster_seqs, threshlod):
	cmdline = [HIERARCH_PATH, "-c", str(threshlod), "-q"]
	child = subprocess.Popen(cmdline, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
	fr.write_fasta(cluster_seqs, child.stdin)
	child.stdin.close()
	clusters = fr.read_cluster(child.stdout)
	return clusters


def split_cluster(cluster_seqs, threshlod):
	if len(cluster_seqs) == 1:
		return {"Cluster" : cluster_seqs}
	
	clusters = []
	preclusters = run_graph(cluster_seqs, threshlod)
	for cl_seqs in preclusters.itervalues():
		hier_clusters = run_hierarchial(cl_seqs, threshlod)
		clusters += hier_clusters.values()
	return {"Clust{0}".format(n) : seqs for n, seqs in enumerate(clusters)}


def correct_reads(cluster_stream, threshlod, out_stream):
	out_seqs = {}
	count = 0
	init_clusters = fr.read_cluster(cluster_stream)
	for cl_name, cl_seqs in init_clusters.iteritems():
		#sys.stderr.write(str(len(cl.seqs)) + " ")
		clusters = split_cluster(cl_seqs, threshlod)
		for newc_name, newc_seqs in clusters.iteritems():
			cons = msa.get_consensus(newc_seqs.keys(), newc_seqs)
			size = 0
			for s in newc_seqs:
				header = s.split(" ")[0]
				size += int(header.split("_")[1]) 
			out_seqs["Seq_{0}_{1}".format(count, size)] = cons
			count += 1
			#sys.stderr.write(str(count) + "\n")
	fr.write_fasta(out_seqs, out_stream)


def main():
	THRESHOLD = 4
	if len(sys.argv) < 2:
		print "USAGE: correct_reads.py reads_clusters"
		return

	correct_reads(open(sys.argv[1], "r"), THRESHOLD, sys.stdout)

if __name__ == "__main__":
	main()
