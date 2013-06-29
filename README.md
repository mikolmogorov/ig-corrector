Llama-fixer
===========

DESCRIPTION
-----------

Llama-fixer is a tool for immunoglobulines NGS data processing.
It provides both duplicates removal and error correction. For quick
manual see "USAGE" section, or see entire file for detailed pipeline
description.


INSTALL
-------

Before building, ensure that you have all dependencies installed.

    make

DEPENDENCIES
------------

* muscle (should be in your executable path)
* Biopython
* Cython
* NumPy

USAGE
-----

    llama-fixer.py [-h] [-q] -o out_dir [-c config] reads_file

    positional arguments:
      reads_file  Input file with reads (fasta or fastq)

    optional arguments:
      -h, --help  show this help message and exit
      -q          Use fastq instead of fasta
      -o out_dir  Output directory
      -c config   Config file (default: config.json)


DATA REQUIREMENTS
-----------------

Reads should be obtained by 454 amplicon sequencing protocol with 2-side
reading. They are expected to fully cover amplicons (immunoglobuline regions).
Also, piplne relies on CDR3 region extracion, so some markers for this
extracion should be provided (see below).


CONFIG DESCRIPTION
------------------

Sample config file (in json format):

    [
        {
            "sampleName" : "VH",
            "forwardMid" : "ATAGATAGAC", 
            "reverseMid" : "ATATAGTCGC",
            "forwardPrimer" : "GCCTACGGCAGCCGCTGGATTGTTATTAC",
            "reversePrimer" : "CACAGACGGGCCTTTTGTAGAC",
            "cdr3Start" : ["YYC"],
            "cdr3End" : ["WG[KQR]"],
            "cdr3Threshold" : 2,
            "sequenceThreshold" : 4 
        },
        {
            "sampleName" : "VL",
            "forwardMid" : "ATCTACTGAC", 
            "reverseMid" : "CACGTAGATC",
            "forwardPrimer" : "GCTGCTGCTGGTCTGCTGCTCCTCGCTG", 
            "reversePrimer" : "GGCGGGAAAATAAAAACAGACGG",
            "cdr3Start" : ["YYC"],
            "cdr3End" : ["GT[KQ]"],
            "cdr3Threshold" : 2,
            "sequenceThreshold" : 4 
        }
    ]

Parameters desctiption:

    sampleName: just name you prefer
    forwardMid: MID sequence for forward read
    reverseMid: sampe for reverse
    forwardPrimer: primer seqence for forward read (see description of "separate.py")
    reversePrimer: same for reverse
    cdr3Start: characteristic protein sequence for start of cdr3 region (see description of "find\_cdr3.py")
    cdr3End: same for end
    cdr3Threshold: threshold for "graph\_clust" for cdr clustering (see pipeline description)
    sequenceThreshold: threshold for "graph\_clust" for reads clustering


PIPELINE FLOW
-------------

0. Separate data by MID, primer clipping ("separate.py")
1. Remove duplicate reads ("remove\_duplicates.py")
2. Extract cdr3 regions ("find\_cdr3.py")
3. Cluster cdr3 ("graph\_clust")
4. Correct cdr3 ("correct\_cdr.py")
5. Correct reads ("correct\_reads.py")


PIPELINE DESCRIPTION
--------------------

### C++ tools

#### graph\_clust

Clustering seqences by edit distance. Two vertexes (seqences) are adjacent, if their
edit distance is less than chosen threshold. Algorithm outputs connected components
in this graph in "cluster" file format (see below). To speed up graph construction, 
preclustering by k-mers is performed (same as BLAST's "anchors"). Choice of k doesn`t 
affect results, but running time. For read length of 450, k = 21 is a good choice.

    Usage: graph_clust -k kmer_size -m num_missmatch [-q] reads_file
    If reads_file is not set, reading from standard input
    -q: quiet

#### hierarchial\_clust

Perform hierarchial sequence clustering. At first, all seqences represents separate
clusters. Then, at each step of algorithm, clusters with the least distance are merged 
(distance between cluster is an average distance between all pairs of seqences).
This continues until distance is less than chosen cutoff. Output is in "cluster"
file format (see below).

    Usage: hierarchial_clust -c cutoff [-q] reads_file
    If reads_file is not set, reading from standard input
    -q: quiet

#### xalign

Protein (as querry) vs nucleotide (as database) alignment. Algorithm is based on
network alignment and supports amino-acid variations. For example: "WG[KQR]" means
W and G on positions 1 and 2, and K, Q or R on position 3. Outputs only alignments
with score >= threshold.

    Usage: xalign [-m match_score] [-x missmatch_score] [-i insert_score]
                  [-d delete_score] -t threshold -q querry reads_file
    Default values: m = 2; x = -2; d = -1; i = -1
    If reads_file is not set, reading from standard input

    Output format: >seq_id (start1, end1, score1) (start2, end2, score2) (...) ...


### Python scripts

#### separate.py

Separate seqencing input data (either fasta or fastq) on samples by MID. Also, perform
primer clipping (if primer sequence provided). Since we have 2-side reading, MIDs
and primers from reverse strand are also clipped. If MID and primer sequence are
not found on forward strand, read goes to "filtered.fasta". If forward MID and 
primer seqence are found on forward strand, but they not persist in reverse strand, 
read is threated as incomplete and goes to "incomplete.fasta" in corresponding 
directory.

#### find\_cdr3.py

Extract cdr3 regions from seqences. Regions are determined by characteristic begin
and end (should be provided as a protein sequence). Algorithm uses "xalign" to find
all possible begins and ends and then performs some magic heuristics.

#### correct\_cdr.py

Correct cdr3 region seqences, that are previously extracted by "find\_cdr3.py". Script
accepts "cluster" file with cdrs (clustering is performed with "graph\_clust"). Similar 
cdrs inside cluster (and seqences, which they correspond to) are merged according to 
next rules: Firstly, indels are corrected (since they are very unlikely to be a real
mutations). Secondly, cdrs with 1-2 missmatches are merdged, if they represent
similar sets of seqences.

#### correct\_reads.py

Perform final read correction. Script accepts corrected cdr "cluster" file. Reads inside
each cluster are being hierarchally clustered (with "hierarchial\_clust"), and then
consensus of each new subcluster is taken. Each subcluster represents a separate real
amplicon and consnsus represents corrected seqence of it (in case of coverage > 1).


FILE TYPES USED
---------------

### Cluster file format (.cl)

Ordinary fasta, separated by cluster names, started with =.

    =Cluster_843_1
    >Seq14996_1
    TATTACTGTGCGAGAGTTGTTACAGGCGGCTCTGAAGACTAC
    =Cluster_844_3
    >Seq25066_1
    TATTACTGTGCGAGAAGGCTCGATACAGACTATGGAAATGACTATGACGTC
    >Seq19110_1
    TATTACTGTGCGAGAAGGCTCGATACAGACTATGGAAATGACTATGACGTC
    >Seq12219_1
    TATTACTGTGCGAGAAGGCTCGATACAGACTATGGAAATGACTATGACGTC


AUTHORS
-------

Kolmogorov Mikhail, Saint-Petersburg Academic University, fenderglass@gmail.com

made for Biocad pharma company.

THIRD PARTY
-----------

* align python module (included): https://github.com/FredrikAppelros/align
* muscle alignment tool
* Biopython
* Cython
* Numpy
