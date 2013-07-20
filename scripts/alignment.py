import sys
import numpy
import logging
import subprocess
from cStringIO import StringIO

import fasta_reader as fr
import align
import ext_tools


logger = logging.getLogger(__name__)


class Aligner:
    def __init__(self, match, missmatch, gap_open, gap_extend):
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.subs = (numpy.ones((256, 256)) * missmatch +
                        numpy.identity(256) * (match - missmatch))
        self.subs = self.subs.astype(numpy.int16)

    def align(self, seq1, seq2, local = False):
        s1 = align.string_to_alignment(seq1)
        s2 = align.string_to_alignment(seq2)
        score, a1, a2 = align.align(s1, s2, self.gap_open, self.gap_extend,
                                    self.subs, local)
        res1, res2 = align.alignment_to_string(a1), align.alignment_to_string(a2)

        if local:
            strip1, strip2 = res1.replace("-", ""), res2.replace("-", "")
            start1, start2 = seq1.index(strip1), seq2.index(strip2)
            start_flank = max(start1, start2)
            end_flank = max(len(seq1) - len(strip1) - start1,
                            len(seq2) - len(strip2) - start2)
            res1 = "-" * start_flank + res1 + "-" * end_flank
            res2 = "-" * start_flank + res2 + "-" * end_flank
        return res1, res2, score

    def edit_dist(self, seq1, seq2):
        al = self.align(seq1, seq2)
        miss = 0
        for a, b in zip(al[0], al[1]):
            if a != b: miss += 1
        return miss

def align_muscle(headers, seqs):
    fasta_dict = {h: seqs[h] for h in headers}
    return ext_tools.muscle(fasta_dict)


def get_consensus(headers, seqs):
    if len(headers) == 1:
        return seqs[headers[0]]

    align = align_muscle(headers, seqs)

    seq_len = len(align[align.keys()[0]])
    result = ""
    for i in xrange(seq_len):
        freq = {}
        for h in align:
            vals = h.split("_")
            mult = int(vals[1]) if len(vals) > 1 else 1
            freq[align[h][i]] = freq.get(align[h][i], 0) + mult

        n = max(freq, key = freq.get)
        if n != "-":
            result += n
    return result
