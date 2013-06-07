#!/usr/bin/env python

from scripts.find_cdr3 import find_cdr3
from scripts.correct_cdr import correct_cdr
from scripts.correct_reads import correct_reads
import sys
import os
import subprocess

BINARIES_PATH = "src"
GRAPH_CLUST_EXEC = "graph_clust"

def cluster_cdr(in_stream, out_stream):
	TRHLD = 2
	cmdline = [GRAPH_CLUST_EXEC, "-k", "11", "-m", str(TRHLD)]
	child = subprocess.Popen(cmdline, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
	for line in in_stream:
		child.stdin.write(line)
	child.stdin.close()
	for line in child.stdout:
		out_stream.write(line)


def process_sample(samplepref, filename, outdir):
	CDR_FILE = outdir + "/" + samplepref + "_cdr.fasta"
	CDR_CLUST_FILE = outdir + "/" + samplepref + "_cdr.cl"
	CDR_CORR_FILE = outdir + "/" + samplepref + "_cdr_corrected.cl"
	READ_CORR_FILE = outdir + "/" + samplepref + "_corrected.fasta"

	#extract cdr3 regions
	CDR3_START = ["YYC"]
	CDR3_END = ["WG[QKR]"] 	#for VH
	#CDR3_END = ["GT[KQ]"]		#for VK, VL
	CDR_TRHLD = 15
	find_cdr3(open(filename, "r"), CDR3_START, CDR3_END, 
				CDR_TRHLD, open(CDR_FILE, "w"))

	#cluster cdr3s
	cluster_cdr(open(CDR_FILE, "r"), open(CDR_CLUST_FILE, "w"))

	#correct cdr
	CDR_CORR_TRHLD = 2
	correct_cdr(open(CDR_CLUST_FILE, "r"), open(filename, "r"), 
				CDR_CORR_TRHLD, open(CDR_CORR_FILE, "w"))

	#correct reads
	READ_CORR_TRHLD = 4
	correct_reads(open(CDR_CORR_FILE, "r"), 
					READ_CORR_TRHLD, open(READ_CORR_FILE, "w"))

def main():
	os.environ["PATH"] += os.pathsep + BINARIES_PATH

	if len(sys.argv) < 4:
		print "USAGE: llama-fixer.py sample_name sample_fasta out_dir"
		return

	process_sample(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
	main()
