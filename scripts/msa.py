import subprocess
import fasta_reader as fr
import logging

MUSCLE_PATH = "muscle"

def align_muscle(headers, seqs):
    fasta_dict = {h: seqs[h] for h in headers}
    logging.getLogger(__name__).debug("Running muscle for {0} seqs"
                                                .format(len(fasta_dict)))

    cmdline = [MUSCLE_PATH, "-diags", "-maxiters", "2", "-quiet"]
    child = subprocess.Popen(cmdline, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    fr.write_fasta(fasta_dict, child.stdin)
    child.stdin.close()
    out_dict = fr.read_fasta(child.stdout)
    logging.getLogger(__name__).debug("Muscle finished")
    return out_dict


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

