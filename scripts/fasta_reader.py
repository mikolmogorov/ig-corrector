def fasta_source(filename):
	fd = open(filename, "r")
	seq = ""
	header = ""
	for line in fd:
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
	fd.close()

def get_seqs(filename):
	seqs = {}
	for h, seq in fasta_source(filename):
		seqs[h] = seq
	return seqs

def iter_kmers(seq, k):
    n = len(seq)
    for i in xrange(0, n - k + 1):
        yield i, seq[i:i + k]
