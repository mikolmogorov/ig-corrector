#!/usr/bin/env python

import sys
import subprocess
import fasta_reader as fr
import msa
import logging


GRAPH_EXEC = "graph_clust"
HIERARCH_EXEC = "hierarchial_clust"


def run_graph(cluster_seqs, threshlod):
    K = 21
    logging.getLogger(__name__).debug("graph_clust started with k = {0} and m = {1}"
                                                            .format(K, threshlod))
    cmdline = [GRAPH_EXEC, "-k", str(K), "-m", str(threshlod), "-q"]
    child = subprocess.Popen(cmdline, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    fr.write_fasta(cluster_seqs, child.stdin)
    child.stdin.close()
    preclusters = fr.read_cluster(child.stdout)
    logging.getLogger(__name__).debug("graph_clust finished")
    return preclusters


def run_hierarchial(cluster_seqs, threshlod):
    logging.getLogger(__name__).debug("hierarchial_clust started with cutoff {0}"
                                                                .format(threshlod))
    cmdline = [HIERARCH_EXEC, "-c", str(threshlod), "-q"]
    child = subprocess.Popen(cmdline, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    fr.write_fasta(cluster_seqs, child.stdin)
    child.stdin.close()
    clusters = fr.read_cluster(child.stdout)
    logging.getLogger(__name__).debug("hierarchial_clust finished")
    return clusters


def split_cluster(cluster_seqs, threshlod):
    logging.getLogger(__name__).debug("Spliiting cluster with size {0}"
                                                        .format(len(cluster_seqs)))
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
    logging.getLogger(__name__).info("Correcting reads with {0} clusters and threshold {1}"
                                                    .format(len(init_clusters), threshlod))
    for cl_name, cl_seqs in init_clusters.iteritems():
        clusters = split_cluster(cl_seqs, threshlod)
        for newc_name, newc_seqs in clusters.iteritems():
            cons = msa.get_consensus(newc_seqs.keys(), newc_seqs)
            size = 0
            for s in newc_seqs:
                header = s.split(" ")[0]
                size += int(header.split("_")[1])
            out_seqs["Seq_{0}_{1}".format(count, size)] = cons
            count += 1
    fr.write_fasta(out_seqs, out_stream)
    logging.getLogger(__name__).info("Correcting finished, {0} unique seqs found"
                                                            .format(len(out_seqs)))


def main():
    logging.basicConfig(level = logging.DEBUG)
    THRESHOLD = 4
    if len(sys.argv) < 2:
        print "USAGE: correct_reads.py reads_clusters"
        return

    correct_reads(open(sys.argv[1], "r"), THRESHOLD, sys.stdout)

if __name__ == "__main__":
    main()
