#!/usr/bin/env python

import sys
import json
import logging
import os
import fasta_reader as fr
import alignment as aln
from Bio.Seq import Seq
from Bio import SeqIO
from collections import namedtuple


PrimerPair = namedtuple("PrimerPair", ["forward", "reverse"])

logger = logging.getLogger(__name__)

MATCH = 1
MISSMATCH = -1
GAP_OPEN = -1
GAP_EXT = -1
aligner = aln.Aligner(MATCH, MISSMATCH, GAP_OPEN, GAP_EXT)

class AlignRes:
    ok = 0
    incomplete = 1
    nomid = 2


def local_alignment(seq1, seq2):
    return aligner.align(seq1, seq2, True)


def with_sequence(mid_pair, primer_pair, seq):
    FWD_MISSMATCH = 6
    REV_MISSMATCH = 10

    rev_seq = str(Seq(seq).reverse_complement())

    for f_dir, op_dir in [(0, 1), (1, 0)]:
        #speeding up alignments
        start_flank = 2 * (len(mid_pair[f_dir]) + len(primer_pair[f_dir]))
        end_flank = 4 * (len(mid_pair[op_dir]) + len(primer_pair[op_dir]))

        start_seq = mid_pair[f_dir] + primer_pair[f_dir]
        start_aln_a, start_aln_b, start_aln_score = local_alignment(seq[0 : start_flank], start_seq)
        if start_aln_score < len(start_seq) - FWD_MISSMATCH:
            continue

        end_seq = mid_pair[op_dir] + primer_pair[op_dir]
        end_aln_a, end_aln_b, end_aln_score = local_alignment(rev_seq[0 : end_flank], end_seq)

        if end_aln_score < len(end_seq) - REV_MISSMATCH:
            return None, AlignRes.incomplete

        start = len(start_aln_b.rstrip("-"))
        end = len(end_aln_b.rstrip("-"))

        res = seq[start : -end]
        if f_dir == 0:
            return res, AlignRes.ok
        else:
            return str(Seq(res).reverse_complement()), AlignRes.ok

    return None, AlignRes.nomid


def get_seqs(filename, is_fastq):
    file = open(filename, "r")
    if is_fastq:
        records = SeqIO.parse(file, "fastq")
    else:
        records = SeqIO.parse(file, "fasta")
    return {rec.id : str(rec.seq) for rec in records}


def split(config_name, reads_file, out_dir, is_fastq):
    logger.info("Splitting {0} started".format(reads_file))
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

    reads = get_seqs(reads_file, is_fastq)
    old_percent = -1
    counter = 0

    for header, seq in reads.iteritems():
        counter += 1
        perc = 10 * counter / len(reads)
        if perc != old_percent:
            old_percent = perc
            sys.stderr.write(str(perc) + " ")

        stripped = seq.strip("acgtn")
        if "N" in stripped:
            logger.debug("{0}: contains N`s".format(header))
            fr.write_fasta({header : seq}, filtered)
            continue

        match_count = 0
        for sample in mids.keys():
            res_seq, res = with_sequence(mids[sample], primers[sample], stripped)
            if res == AlignRes.ok:
                logger.debug("{0}: is {1}".format(header, sample))
                fr.write_fasta({header : res_seq}, files[sample])
                match_count += 1
                break
            elif res == AlignRes.incomplete:
                logger.debug("{0}: is incomplete {1}".format(header, sample))
                fr.write_fasta({header : stripped}, incomplete[sample])
                match_count += 1
                break

        if match_count == 0:
            logger.debug("{0}: no mid found".format(header))
            fr.write_fasta({header : seq}, filtered)

    sys.stderr.write("\n")
    logger.info("Splitting finished")


def main():
    #logging.basicConfig(level = logging.DEBUG)
    split(sys.argv[1], sys.argv[2], sys.argv[3], True)


if __name__ == "__main__":
    main()
