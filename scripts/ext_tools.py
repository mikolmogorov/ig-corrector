import logging
import subprocess
import fasta_reader as fr
from cStringIO import StringIO

logger = logging.getLogger(__name__)

GRAPH_CLUST_EXEC = "graph_clust"
HIERARCH_CLUST_EXEC = "hierarchial_clust"
MUSCLE_EXEC = "muscle"

def graph_clust(fasta_dict, kmer, missmatches):
    logger.debug(GRAPH_CLUST_EXEC + " for {0} seqs started".format(len(fasta_dict)))

    buffer = StringIO()
    fr.write_fasta(fasta_dict, buffer)

    cmdline = [GRAPH_CLUST_EXEC, "-k", str(kmer), "-m", str(missmatches), "-q"]
    process = subprocess.Popen(cmdline, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    child_stdout, _ = process.communicate(input=buffer.getvalue())

    logger.debug(GRAPH_CLUST_EXEC + " finished")
    return fr.read_cluster(StringIO(child_stdout))


def hierarchial_clust(fasta_dict, cutoff):
    logger.debug(HIERARCH_CLUST_EXEC + " for {0} seqs started".format(len(fasta_dict)))

    buffer = StringIO()
    fr.write_fasta(fasta_dict, buffer)

    cmdline = [HIERARCH_CLUST_EXEC, "-c", str(cutoff), "-q"]
    proc = subprocess.Popen(cmdline, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    child_stdout, _ = proc.communicate(input=buffer.getvalue())

    logger.debug(HIERARCH_CLUST_EXEC + " finished")
    return fr.read_cluster(StringIO(child_stdout))


def muscle(fasta_dict):
    logger.debug("Running muscle for {0} seqs".format(len(fasta_dict)))

    buffer = StringIO()
    fr.write_fasta(fasta_dict, buffer)

    cmdline = [MUSCLE_EXEC, "-diags", "-maxiters", "2", "-quiet"]
    proc = subprocess.Popen(cmdline, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    child_stdout, _ = proc.communicate(input=buffer.getvalue())

    logger.debug("Muscle finished")
    return fr.read_fasta(StringIO(child_stdout))

