import sys, os
import numpy

sys.path.insert(1, os.path.join(os.path.dirname(__file__),"../third-party/calign"))
import align


class Aligner:
    def __init__(self, match, missmatch, gap_open, gap_extend):
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.subs = (numpy.ones((256, 256)) * missmatch +
                        numpy.identity(256) * (match - missmatch))
        self.subs = self.subs.astype(numpy.int16)

    def align(self, seq1, seq2, local = False):
        s1 = calign.string_to_alignment(seq1)
        s2 = calign.string_to_alignment(seq2)
        score, a1, a2 = calign.align(s1, s2, self.gap_open, self.gap_extend,
                                    self.subs, local)
        res1, res2 = calign.alignment_to_string(a1), calign.alignment_to_string(a2)

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


def main():
    a = Aligner(0, -1, -1, -1)
    print a.align(sys.argv[1], sys.argv[2], False)


if __name__ == "__main__":
    main()
