#!/usr/bin/env python

import sys, os
import fasta_reader as fr
import alignment
import logging
import ext_tools


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


KMER = 11
kmer_cache = KmerCache(KMER)


def run_graph(cluster_seqs, threshold):
    GRAPH_KMER = 21
    return ext_tools.graph_clust(cluster_seqs, GRAPH_KMER, threshold)


def run_hierarchial(cluster_seqs, threshold):
    return ext_tools.hierarchial_clust(cluster_seqs, threshold)


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
        _, qty = header.split("_")
        real_seqs += int(qty)

    #naive kmer correction for 2-sequence consensus
    if real_seqs > 2:
        cons = alignment.get_consensus(seqs.keys(), seqs)
        kmer_cache.add(cons)
    else:
        aln = alignment.align_muscle(seqs.keys(), seqs)
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


def correct_reads(cluster_stream, threshlod, out_corr_stream, out_cluster_stream):
    out_seqs = {}
    count = 0
    init_clusters = fr.read_cluster(cluster_stream)
    logger.info("Correcting reads with {0} clusters and threshold {1}"
                                        .format(len(init_clusters), threshlod))
    for cl_name, cl_seqs in init_clusters.iteritems():
        clusters = split_cluster(cl_seqs, threshlod)
        fr.write_cluster(clusters, out_cluster_stream)
        for newc_name, newc_seqs in clusters.iteritems():
            cons_hdr, cons_seq = get_consensus(newc_seqs, count)
            out_seqs[cons_hdr] = cons_seq
            count += 1
    fr.write_fasta(out_seqs, out_corr_stream)
    logger.info("Correcting finished, {0} unique seqs found".format(len(out_seqs)))


def main():
    logging.basicConfig(level = logging.DEBUG)
    THRESHOLD = 4
    if len(sys.argv) < 2:
        print "USAGE: correct_reads.py reads_clusters"
        return

    correct_reads(open(sys.argv[1], "r"), THRESHOLD, sys.stdout, open(os.devnull, "w"))

if __name__ == "__main__":
    main()
