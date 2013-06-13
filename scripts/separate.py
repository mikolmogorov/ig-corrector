#!/usr/bin/env python2

import sys
import json
import logging
import os
import fasta_reader as fr
from Bio.Seq import Seq
from Bio import pairwise2
from collections import namedtuple


PrimerPair = namedtuple("PrimerPair", ["forward", "reverse"])


class AlignRes:
    ok = 0
    incomplete = 1
    nomid = 2


def with_sequence(mid_pair, primer_pair, seq):
    FWD_MISSMATCH = 6
    REV_MISSMATCH = 10

    rev_seq = str(Seq(seq).reverse_complement())

    for f_dir, op_dir in [(0, 1), (1, 0)]:
        start_flank = 2 * (len(mid_pair[f_dir]) + len(primer_pair[f_dir]))
        end_flank = 4 * (len(mid_pair[op_dir]) + len(primer_pair[op_dir]))

        start_seq = mid_pair[f_dir] + primer_pair[f_dir]
        start_aln = pairwise2.align.localms(seq[0 : start_flank], start_seq,
                                            1, -1, -1, -1, one_alignment_only = 1)[0]
        if start_aln[2] < len(start_seq) - FWD_MISSMATCH:
            continue

        end_seq = mid_pair[op_dir] + primer_pair[op_dir]
        end_aln = pairwise2.align.localms(rev_seq[0 : end_flank], end_seq, 1, -1, -1, -1,
                                            one_alignment_only = 1)[0]
        if end_aln[2] < len(end_seq) - REV_MISSMATCH:
            return None, AlignRes.incomplete

        start = len(start_aln[1].rstrip("-"))
        end = len(end_aln[1].rstrip("-"))

        res = seq[start : -end]
        if f_dir == 0:
            return res, AlignRes.ok
        else:
            return str(Seq(res).reverse_complement()), AlignRes.ok

    return None, AlignRes.nomid


def split(config_name, reads_file, out_dir):
    logging.getLogger(__name__).info("Splitting {0} started".format(reads_file))
    config = json.load(open(config_name, "r"))
    mids = {}
    primers = {}
    files = {}
    incomplete = {}
    filtered = open(os.path.join(out_dir, "filtered.fasta"), "w")
    for sample in config:
        name = sample[u"sampleName"].encode("ascii")
        sample_dir = os.path.join(out_dir, name)
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)

        fwdMid = sample[u"forwardMid"].encode("ascii")
        revMid = sample[u"reverseMid"].encode("ascii")
        fwdPrimer = sample[u"forwardPrimer"].encode("ascii")
        revPrimer = sample[u"reversePrimer"].encode("ascii")
        mids[name] = PrimerPair(fwdMid, revMid)
        primers[name] = PrimerPair(fwdPrimer, revPrimer)
        files[name] = open(os.path.join(sample_dir, name + ".fasta"), "w")
        incomplete[name] = open(os.path.join(sample_dir, name + "_incomplete.fasta"), "w")

    processed = 0
    for header, seq, _ in fr._fastq_source(open(reads_file, "r")):
        print processed
        processed += 1

        stripped = seq.strip("acgtn")
        if "N" in stripped:
            logging.getLogger(__name__).debug("{0}: contains N`s".format(header))
            fr.write_fasta({header : seq}, filtered)
            continue

        match_count = 0
        for sample in mids.keys():
            res_seq, res = with_sequence(mids[sample], primers[sample], stripped)
            if res == AlignRes.ok:
                logging.getLogger(__name__).debug("{0}: is {1}".format(header, sample))
                fr.write_fasta({header : res_seq}, files[sample])
                match_count += 1
                break
            elif res == AlignRes.incomplete:
                logging.getLogger(__name__).debug("{0}: is incomplete {1}".format(header, sample))
                fr.write_fasta({header : stripped}, incomplete[sample])
                match_count += 1
                break

        if match_count == 0:
            logging.getLogger(__name__).debug("{0}: no mid found".format(header))
            fr.write_fasta({header : seq}, filtered)

    logging.getLogger(__name__).info("Splitting finished")


def main():
    #logging.basicConfig(level = logging.DEBUG)
    split(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()
