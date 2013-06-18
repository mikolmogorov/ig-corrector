#!/usr/bin/env python

import fasta_reader as fr
import sys
import logging

logger = logging.getLogger(__name__)

def remove_dups(stream_in, stream_out):
    seqs = fr.read_fasta(stream_in)
    logger.info("Removing duplicates from {0} sequences".format(len(seqs)))
    counter = 0
    new_seqs = {}
    seq_count = {}
    for h in seqs:
        seq_count[seqs[h]] = seq_count.get(seqs[h], 0) + 1
    for s in seq_count:
        new_seqs["Seq{0}_{1}".format(counter, seq_count[s])] = s
        counter += 1
    fr.write_fasta(new_seqs, stream_out)
    logger.info("Done, {0} sequences left".format(len(new_seqs)))

def main():
    if len(sys.argv) < 2:
        print "USAGE: remove_duplicates.py filename"
        return

    newSeqs = remove_dups(open(sys.argv[1], "r"), sys.stdout)

if __name__ == "__main__":
    main()
