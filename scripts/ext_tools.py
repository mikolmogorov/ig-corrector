import logging
import subprocess
import tempfile
import fasta_reader as fr
from collections import namedtuple
import re


GRAPH_CLUST_EXEC = "graph_clust"
HIERARCH_CLUST_EXEC = "hierarchial_clust"
MUSCLE_EXEC = "muscle"
XALIGN_EXEC = "xalign"


logger = logging.getLogger(__name__)
AlignInfo = namedtuple("AlignInfo", ["begin", "end", "score"])


def graph_clust(fasta_dict, kmer, missmatches):
    logger.debug(GRAPH_CLUST_EXEC + " for {0} seqs started".format(len(fasta_dict)))

    in_file = tempfile.SpooledTemporaryFile()
    out_file = tempfile.SpooledTemporaryFile()
    fr.write_fasta(fasta_dict, in_file)
    in_file.flush()
    in_file.seek(0)

    cmdline = [GRAPH_CLUST_EXEC, "-k", str(kmer), "-m", str(missmatches), "-q"]
    process = subprocess.Popen(cmdline, stdin=in_file, stdout=out_file)
    process.wait()
    out_file.flush()
    out_file.seek(0)

    logger.debug(GRAPH_CLUST_EXEC + " finished")
    return fr.read_cluster(out_file)


def hierarchial_clust(fasta_dict, cutoff):
    logger.debug(HIERARCH_CLUST_EXEC + " for {0} seqs started".format(len(fasta_dict)))

    in_file = tempfile.SpooledTemporaryFile()
    out_file = tempfile.SpooledTemporaryFile()
    fr.write_fasta(fasta_dict, in_file)
    in_file.flush()
    in_file.seek(0)

    cmdline = [HIERARCH_CLUST_EXEC, "-c", str(cutoff), "-q"]
    proc = subprocess.Popen(cmdline, stdin=in_file, stdout=out_file)
    proc.wait()
    out_file.flush()
    out_file.seek(0)

    logger.debug(HIERARCH_CLUST_EXEC + " finished")
    return fr.read_cluster(out_file)


def xalign(fasta_dict, query, match, missmatch, indel):
    CODON_LEN = 3
    qry_len = len(re.sub(r"\[.*\]", "X", query)) * CODON_LEN
    trhld = qry_len * match - (match + indel)                   #one indel is allowed
    #assert qry_len == 9

    in_file = tempfile.SpooledTemporaryFile()
    out_file = tempfile.SpooledTemporaryFile()
    fr.write_fasta(fasta_dict, in_file)
    in_file.flush()
    in_file.seek(0)

    logger.debug("Calling xalign with querry \"{0}\"".format(query))
    cmdline = [XALIGN_EXEC, "-t", str(trhld), "-q", query, "-m", str(match),
                    "-x", str(missmatch), "-i", str(indel), "-d", str(indel)]
    child = subprocess.Popen(cmdline, stdin=in_file, stdout=out_file)
    child.wait()
    out_file.flush()
    out_file.seek(0)

    out_dict = {}
    for line in out_file:
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


def muscle(fasta_dict):
    logger.debug("Running muscle for {0} seqs".format(len(fasta_dict)))

    in_file = tempfile.SpooledTemporaryFile()
    out_file = tempfile.SpooledTemporaryFile()
    fr.write_fasta(fasta_dict, in_file)
    in_file.flush()
    in_file.seek(0)

    cmdline = [MUSCLE_EXEC, "-diags", "-maxiters", "2", "-quiet"]
    proc = subprocess.Popen(cmdline, stdin=in_file, stdout=out_file)
    proc.wait()
    out_file.flush()
    out_file.seek(0)

    logger.debug("Muscle finished")
    return fr.read_fasta(out_file)

