#!/usr/bin/env python

import fasta_reader as fr
import msa
import sys
import editdist
import logging
from Bio import pairwise2
from itertools import product, combinations
from collections import defaultdict


logger = logging.getLogger(__name__)


def global_alignment(seq1, seq2):
    MATCH = 0
    MISSMATCH = -1
    GAP_OPEN = -2
    GAP_EXT = -2
    align = pairwise2.align.globalms(seq1, seq2, MATCH, MISSMATCH, GAP_OPEN, GAP_EXT,
                                    one_alignment_only = True)[0]
    return align[0:2]


def choose_true(cdr1, cdr2, weight, threshold):
    if float(weight[cdr1]) / weight[cdr2] >= threshold:
        true_cdr, false_cdr = cdr1, cdr2
    elif len(cdr1) % 3 == 0 and len(cdr2) % 3 != 0:
        true_cdr, false_cdr = cdr1, cdr2
    elif len(cdr2) % 3 == 0 and len(cdr1) % 3 != 0:
        true_cdr, false_cdr = cdr2, cdr1
    else:
        logger.debug("Ambigious cdrs: {0} and {1}".format(cdr1, cdr2))
        true_cdr, false_cdr = cdr1, cdr2
    logger.debug("{0} (as true) merged with {1}".format(cdr1, cdr2))
    return true_cdr, false_cdr


def correct_indels(cdr_map, weight, full_seqs, threshold):
    logger.debug("Correcting indels")
    sorted_r = sorted(weight, key = weight.get, reverse = True)
    for cdr1, cdr2 in combinations(sorted_r, 2):
        if cdr1 not in cdr_map or cdr2 not in cdr_map:
            continue

        align_a, align_b = global_alignment(cdr1, cdr2)
        n_miss = 0
        for i in xrange(len(align_a)):
            if align_a[i] != align_b[i] and align_a[i] != "-" and align_b[i] != "-":
                n_miss += 1

        if n_miss > 0:
            continue

        true_cdr, false_cdr = choose_true(cdr1, cdr2, weight, threshold)

        cdr_map[true_cdr].extend(cdr_map[false_cdr])
        del cdr_map[false_cdr]
        weight[true_cdr] += weight[false_cdr]
        del weight[false_cdr]


#should be called _after_ indels correction
def correct_missmatches(cdr_map, weight, full_seqs, threshold):
    logger.debug("Correcting missmatches")
    cons_cache = {}
    sorted_r = sorted(weight, key = weight.get, reverse = True)
    for cdr1, cdr2 in combinations(sorted_r, 2):
        if cdr1 not in cdr_map or cdr2 not in cdr_map:
            continue

        if cdr1 not in cons_cache:
            cons_cache[cdr1] = msa.get_consensus(cdr_map[cdr1], full_seqs)
        if cdr2 not in cons_cache:
            cons_cache[cdr2] = msa.get_consensus(cdr_map[cdr2], full_seqs)
        dist = editdist.distance(cons_cache[cdr1], cons_cache[cdr2])
        if dist > 4:
            continue

        true_cdr, false_cdr = choose_true(cdr1, cdr2, weight, threshold)

        cdr_map[true_cdr].extend(cdr_map[false_cdr])
        del cdr_map[false_cdr]
        weight[true_cdr] += weight[false_cdr]
        del weight[false_cdr]

        cons_cache.pop(true_cdr, None)
        cons_cache.pop(false_cdr, None)


def correct_cluster(clust_seqs, full_seqs, threshold):
    logger.debug("Correcting cluster of size {0} with threshold {1}"
                                    .format(len(clust_seqs), threshold))
    cdr_map = defaultdict(list)
    weight = {}
    for seq_head, seq_cdr in clust_seqs.iteritems():
        cdr_map[seq_cdr].append(seq_head)
        qty = int(seq_head.split("_")[1])
        weight[seq_cdr] = weight.get(seq_cdr, 0) + qty

    if len(cdr_map) == 1:
        return {"Cluster" : clust_seqs}

    logger.debug("Uniqe cdrs: {0}".format(len(cdr_map)))
    correct_indels(cdr_map, weight, full_seqs, threshold)
    correct_missmatches(cdr_map, weight, full_seqs, threshold)
    logger.debug("Uniqe cdrs after correction: {0}".format(len(cdr_map)))

    clusters = {}
    counter = 0
    for cdr in cdr_map:
        clust_seqs = {}
        for head in cdr_map[cdr]:
            clust_seqs[head] = cdr
        clusters["Cluster{0}".format(counter)] = clust_seqs
        counter += 1
    return clusters


def out_cluster(name, clust_seqs, full_seqs, out_stream):
    out_stream.write("=" + name + "\n")
    for hdr, seq in clust_seqs.iteritems():
        out_stream.write(">{0} {1}\n{2}\n".format(hdr, seq, full_seqs[hdr]))


def correct_cdr(cdr_stream, seqs_stream, threshold, out_stream):
    full_seqs = fr.read_fasta(seqs_stream)
    clusters = fr.read_cluster(cdr_stream)
    logger.info("Cdr correction started, {0} cdr clusters".format(len(clusters)))
    cl_id = 0
    for c_name, c_seqs in clusters.iteritems():
        for newc_name, newc_seqs in correct_cluster(c_seqs, full_seqs, threshold).iteritems():
            out_cluster("Cluster_{0}".format(cl_id), newc_seqs, full_seqs, out_stream)
            cl_id += 1
    logger.info("Cdr correction finished, {0} cdr clusters".format(cl_id))


def main():
    logging.basicConfig(level = logging.DEBUG)
    THRESHOLD = 2
    if len(sys.argv) < 3:
        print "USAGE: correct_cdr.py cdr_cluster_file reads_file"
        return
    correct_cdr(open(sys.argv[1], "r"), open(sys.argv[2], "r"), THRESHOLD, sys.stdout)


if __name__ == "__main__":
    main()
