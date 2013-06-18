#!/usr/bin/env python

import sys
import subprocess
import fasta_reader as fr
import msa
import logging


GRAPH_EXEC = "graph_clust"
HIERARCH_EXEC = "hierarchial_clust"


logger = logging.getLogger(__name__)


class KmerCache:
    def __init__(self, k):
        self.cache = set()
        self.k = k

    def add(self, sequence):
        for i in xrange(len(sequence) - self.k + 1):
            self.cache.add(sequence[i : i + self.k])

    def check(self, kmer):
        return kmer in self.cache

    def count_trusted(self, sequence):
        counter = 0
        for i in xrange(len(sequence) - self.k + 1):
            #print sequence[i : i + self.k]
            if sequence[i : i + self.k] in self.cache:
                counter += 1
        return counter


kmer_cache = KmerCache(11)


def run_graph(cluster_seqs, threshlod):
    K = 21
    logger.debug("graph_clust started with k = {0} and m = {1}".format(K, threshlod))
    cmdline = [GRAPH_EXEC, "-k", str(K), "-m", str(threshlod), "-q"]
    child = subprocess.Popen(cmdline, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    fr.write_fasta(cluster_seqs, child.stdin)
    child.stdin.close()
    preclusters = fr.read_cluster(child.stdout)
    logger.debug("graph_clust finished")
    return preclusters


def run_hierarchial(cluster_seqs, threshlod):
    logger.debug("hierarchial_clust started with cutoff {0}".format(threshlod))
    cmdline = [HIERARCH_EXEC, "-c", str(threshlod), "-q"]
    child = subprocess.Popen(cmdline, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    fr.write_fasta(cluster_seqs, child.stdin)
    child.stdin.close()
    clusters = fr.read_cluster(child.stdout)
    logger.debug("hierarchial_clust finished")
    return clusters


def split_cluster(cluster_seqs, threshlod):
    logger.debug("Spliiting cluster with size {0}".format(len(cluster_seqs)))
    if len(cluster_seqs) == 1:
        return {"Cluster" : cluster_seqs}

    clusters = []
    preclusters = run_graph(cluster_seqs, threshlod)
    for cl_seqs in preclusters.itervalues():
        hier_clusters = run_hierarchial(cl_seqs, threshlod)
        clusters.extend(hier_clusters.values())
    return {"Clust{0}".format(n) : seqs for n, seqs in enumerate(clusters)}


def get_consensus(seqs, seq_id):
    if len(seqs) == 1:
        return "Seq_{0}_1".format(seq_id), seqs.values()[0]

    real_seqs = 0
    for s in seqs:
        header = s.split(" ")[0]
        real_seqs += int(header.split("_")[1])

    #naive kmer correction for 2-sequence consensus
    if real_seqs > 2:
        cons = msa.get_consensus(seqs.keys(), seqs)
        kmer_cache.add(cons)
    else:
        aln = msa.align_muscle(seqs.keys(), seqs)
        seq_a = aln.values()[0]
        seq_b = aln.values()[1]
        cons = ""
        for i in xrange(len(seq_a)):
            if seq_a[i] == seq_b[i]:
                cons += seq_a[i]
            else:
                k = kmer_cache.k
                tr_a = kmer_cache.count_trusted(seq_a[i - k : i + k])
                tr_b = kmer_cache.count_trusted(seq_b[i - k : i + k])
                cons += seq_a[i] if tr_a >= tr_b else seq_b[i]
    return "Seq_{0}_{1}".format(seq_id, real_seqs), cons


def correct_reads(cluster_stream, threshlod, out_stream):
    out_seqs = {}
    count = 0
    init_clusters = fr.read_cluster(cluster_stream)
    logger.info("Correcting reads with {0} clusters and threshold {1}"
                                        .format(len(init_clusters), threshlod))
    for cl_name, cl_seqs in init_clusters.iteritems():
        clusters = split_cluster(cl_seqs, threshlod)
        for newc_name, newc_seqs in clusters.iteritems():
            cons_hdr, cons_seq = get_consensus(newc_seqs, count)
            out_seqs[cons_hdr] = cons_seq
            count += 1
    fr.write_fasta(out_seqs, out_stream)
    logger.info("Correcting finished, {0} unique seqs found".format(len(out_seqs)))


def main():
    logging.basicConfig(level = logging.DEBUG)
    THRESHOLD = 4
    if len(sys.argv) < 2:
        print "USAGE: correct_reads.py reads_clusters"
        return

    correct_reads(open(sys.argv[1], "r"), THRESHOLD, sys.stdout)

if __name__ == "__main__":
    main()
