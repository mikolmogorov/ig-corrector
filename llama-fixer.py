#!/usr/bin/env python

import sys, os
import signal
import subprocess
import logging
import threading
import json
import argparse
import StringIO

BINARIES_PATH = "src"
GRAPH_CLUST_EXEC = "graph_clust"
PYTHON_LIBS = "third-party/calign"

sys.path.append(os.path.abspath(PYTHON_LIBS))
os.environ["PATH"] += os.pathsep + os.path.abspath(BINARIES_PATH)

from scripts.find_cdr3 import find_cdr3
from scripts.correct_cdr import correct_cdr
from scripts.correct_reads import correct_reads
from scripts.separate import split
from scripts.remove_duplicates import remove_dups


log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: %(message)s",
                                                                            "%H:%M:%S")
logger = logging.getLogger()

class LogDistributorHandler(logging.Handler):
    def __init__(self):
        super(LogDistributorHandler, self).__init__()
        self.file_handlers = {}

    def addFileHandler(self, name, file):
        self.file_handlers[name] = file

    def handle(self, record):
        if record.threadName in self.file_handlers:
            self.file_handlers[record.threadName].emit(record)


class SampleThread(threading.Thread):
    def __init__(self, name, *args):
        super(SampleThread, self).__init__()
        self.name = name
        self.args = args

    def run(self):
        process_sample(self.name, *self.args)


class InterruptWatcher:
    def __init__(self):
        self.child = os.fork()
        if self.child == 0:
            return
        else:
            self.watch()

    def watch(self):
        try:
            os.wait()
        except KeyboardInterrupt:
            sys.stderr.write("Keyboard interrupt recived, exiting\n")
            self.kill()
        sys.exit()

    def kill(self):
        try:
            os.kill(self.child, signal.SIGKILL)
        except OSError: pass


def cluster_cdr(in_stream, out_stream, threshold):
    K = 11
    logger.info("graph_clust started with k = {0} and m = {1}".format(K, threshold))

    cmdline = [GRAPH_CLUST_EXEC, "-k", str(K), "-m", str(threshold), "-q"]
    child = subprocess.Popen(cmdline, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    buffer = StringIO.StringIO()
    for line in in_stream:
        buffer.write(line)

    child_stdout, _ = child.communicate(input=buffer.getvalue())
    for line in iter(child_stdout.splitlines()):
        out_stream.write(line + "\n")
    logger.info("graph_clust finished")


def process_sample(samplepref, filename, outdir, cdr_start,
                                cdr_end, cdr_threshold, seq_threshold):
    UNIQUE_FILE = os.path.join(outdir, samplepref + "_unique.fasta")
    CDR_FILE = os.path.join(outdir, samplepref + "_cdr.fasta")
    CDR_CLUST_FILE = os.path.join(outdir, samplepref + "_cdr.cl")
    CDR_CORR_FILE = os.path.join(outdir, samplepref + "_cdr_corrected.cl")
    READ_CORR_FILE = os.path.join(outdir, samplepref + "_corrected.fasta")

    logger.info("Removing duplicate sequqnces...")
    remove_dups(open(filename, "r"), open(UNIQUE_FILE, "w"))

    logger.info("Finding regions...")
    find_cdr3(open(UNIQUE_FILE, "r"), cdr_start, cdr_end, open(CDR_FILE, "w"))

    logger.info("Clustering extracted cdr`s...")
    cluster_cdr(open(CDR_FILE, "r"), open(CDR_CLUST_FILE, "w"), cdr_threshold)

    logger.info("Correcting cdrs...")
    correct_cdr(open(CDR_CLUST_FILE, "r"), open(UNIQUE_FILE, "r"),
                cdr_threshold, open(CDR_CORR_FILE, "w"))

    logger.info("Correcting reads...")
    correct_reads(open(CDR_CORR_FILE, "r"),
                    seq_threshold, open(READ_CORR_FILE, "w"))
    logger.info("Finished!")


def run_jobs(config_name, out_dir, reads_file, log_distr, is_fastq):

    split(config_name, reads_file, out_dir, is_fastq)
    config = json.load(open(config_name, "r"))
    for sample in config:
        name = sample[u"sampleName"].encode("ascii")
        sample_dir = os.path.join(out_dir, name)
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)

        sample_input = os.path.join(sample_dir, name + ".fasta")

        cdr3_start = [e.encode("ascii") for e in sample[u"cdr3Start"]]
        cdr3_end = [e.encode("ascii") for e in sample[u"cdr3End"]]
        cdr3_threshold = int(sample[u"cdr3Threshold"])
        seq_threshold = int(sample[u"sequenceThreshold"])

        log_file = os.path.join(sample_dir, name + "_log.txt")
        file_handler = logging.FileHandler(log_file, mode = "w")
        file_handler.setFormatter(log_formatter)
        log_distr.addFileHandler(name, file_handler)

        thread = SampleThread(name, sample_input, sample_dir, cdr3_start,
                                cdr3_end, cdr3_threshold, seq_threshold)
        thread.start()


def main():
    parser = argparse.ArgumentParser(description = "A data processing utility for immunoglobulines NGS")
    parser.add_argument("reads_file", action = "store", help = "Input file with reads")
    parser.add_argument("-q", action = "store_const", metavar = "fastq", dest = "fastq",
                        default = False, const  = True, help = "Use fastq instead of fasta")
    parser.add_argument("-o", action = "store", metavar = "out_dir", dest = "out_dir",
                        required = True, help = "Output directory")
    parser.add_argument("-c", action = "store", metavar = "config", dest = "config",
                        default = "config.json", help = "Config file (default: config.json)")
    args = parser.parse_args()

    config_file = args.config
    reads_file = args.reads_file
    out_dir = args.out_dir
    is_fastq = args.fastq

    logging.getLogger().setLevel(logging.DEBUG)
    console_formatter = logging.Formatter("[%(asctime)s] %(threadName)s: %(message)s",
                                                                            "%H:%M:%S")
    common_file = os.path.join(out_dir, "main.log")
    common_handler = logging.FileHandler(common_file, mode = "w")
    common_handler.setFormatter(log_formatter)

    console_log = logging.StreamHandler()
    console_log.setLevel(logging.INFO)
    console_log.setFormatter(console_formatter)

    file_distr = LogDistributorHandler()
    file_distr.setLevel(logging.DEBUG)
    file_distr.addFileHandler("MainThread", common_handler)

    logging.getLogger().addHandler(console_log)
    logging.getLogger().addHandler(file_distr)

    InterruptWatcher()

    run_jobs(config_file, out_dir, reads_file, file_distr, is_fastq)


if __name__ == "__main__":
    main()
