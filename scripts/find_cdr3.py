#!/usr/bin/env python

import fasta_reader as fr
import sys
import subprocess
import logging
import re
from collections import namedtuple, defaultdict

XALIGN_EXEC = "xalign"

AlignInfo = namedtuple("AlignInfo", ["begin", "end", "score"])
logger = logging.getLogger(__name__)


def xalign(fasta_dict, query):
    MATCH = 2
    MISSMATCH = -2
    INDEL = -1
    CODON_LEN = 3

    qry_len = len(re.sub(r"\[.*\]", "X", query)) * CODON_LEN
    trhld = qry_len * MATCH - (MATCH + INDEL)       #one indel is allowed
    #assert qry_len == 9

    logger.debug("Calling xalign with querry \"{0}\"".format(query))
    cmdline = [XALIGN_EXEC, "-t", str(trhld), "-q", query, "-m", str(MATCH),
                    "-x", str(MISSMATCH), "-i", str(INDEL), "-d", str(INDEL)]
    child = subprocess.Popen(cmdline, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    fr.write_fasta(fasta_dict, child.stdin)
    child.stdin.close()

    out_dict = {}
    for line in child.stdout:
        values = line.strip().split(" ")
        header = values[0][1:]
        aligns = values[1:]
        #assert header.startswith(">")
        out_dict[header] = []
        for align in aligns:
            start, end, score = align[1:-1].split(",")
            out_dict[header].append(AlignInfo(int(start), int(end), int(score)))
    logger.debug("xalign finished")
    return out_dict


def find_cdr3(in_stream, start_seqs, end_seqs, out_stream):
    MIN_CDR_LEN = 30
    MAX_CDR_LEN = 90

    logger.info("Finding cdr regions started".format(MIN_CDR_LEN, MAX_CDR_LEN))
    seqs = fr.read_fasta(in_stream)
    start_align = defaultdict(list)
    end_align = defaultdict(list)
    for qry in start_seqs:
        for h, alns in xalign(seqs, qry).iteritems():
            start_align[h].extend(alns)
    for qry in end_seqs:
        for h, alns in xalign(seqs, qry).iteritems():
            end_align[h].extend(alns)

    for h, seq in seqs.iteritems():
        candidates = []
        for start in start_align[h]:
            for end in end_align[h]:
                dist = end.begin - start.begin
                if (MIN_CDR_LEN <= dist and dist <= MAX_CDR_LEN and
                                        start.begin > len(seq) / 2):
                    candidates.append(AlignInfo(start.begin, end.begin - 1, start.score + end.score))

        if len(candidates) == 0:
            logger.debug("{0} : no cdr found".format(h))
            continue

        max_score = 0
        cand = None
        for c in candidates:
            if c.score > max_score:
                cand = c
                max_score = c.score

        cdr = seq[cand.begin : cand.end + 1]
        out_stream.write(">{0}\n{1}\n".format(h, cdr))
        logger.debug("{0} : cdr found with score {1}".format(h, max_score))
    logger.info("Finding cdr regions finished")


def main():
    logging.basicConfig(level = logging.DEBUG)
    CDR3_START = "YYC"
    CDR3_END = "WG[QKR]"
    #CDR3_END = "GTK"
    if len(sys.argv) < 2:
        print "USAGE: find_cdr3.py reads_file"
        return

    find_cdr3(open(sys.argv[1], "r"), [CDR3_START], [CDR3_END], sys.stdout)

if __name__ == "__main__":
    main()
