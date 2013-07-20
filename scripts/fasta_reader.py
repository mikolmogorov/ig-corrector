def _fasta_source(stream):
    seq = ""
    header = ""
    for line in stream:
        l = line.strip("\n")
        if l.startswith(">"):
            if len(header) > 0:
                yield header, seq
                seq = ""
            header = l[1:]
        else:
            seq += l
    if seq != "":
        yield header, seq


def _fastq_source(stream):
    seq = ""
    header = ""
    state = 0
    for line in stream:
        l = line.strip("\n")
        if state == 0:
            header = l[1:]
        elif state == 1:
            seq = l
        elif state == 3:
            yield header, seq, line
        state = (state + 1) % 4


def read_fasta(stream):
    seqs = {}
    for h, seq in _fasta_source(stream):
        seqs[h] = seq
    return seqs


def read_fastq(stream):
    seqs = {}
    for h, seq, _ in _fastq_source(stream):
        seqs[h] = seq
    return seqs


def write_fasta(fasta_dict, stream):
    for h in fasta_dict:
        stream.write(">{0}\n{1}\n".format(h, fasta_dict[h]))


def read_cluster(stream):
    clusters = {}

    fasta_buffer = ""
    clust_header = None
    for line in stream:
        if line.startswith("="):
            if clust_header:
                clusters[clust_header] = read_fasta(iter(fasta_buffer.splitlines()))
            clust_header = line.strip()[1:]
            fasta_buffer = ""
        else:
            fasta_buffer += line
    if clust_header:
        clusters[clust_header] = read_fasta(iter(fasta_buffer.splitlines()))
    return clusters


def write_cluster(clusters, stream):
    for cl in clusters:
        stream.write("=" + cl + "\n")
        write_fasta(clusters[cl], stream)
