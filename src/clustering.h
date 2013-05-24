#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <list>
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include "fasta.h"
#include "disjoint_set.h"

typedef std::vector<std::list<FastaRecord>> FastaSet;

class Clusterisator
{
public:
	Clusterisator() {}
	void doJob (FastaSequences& seqs, FastaSet& output, int kmerSize, int nMissmatches);

private:
	typedef std::unordered_map<int, std::string> FastaHash;
	typedef std::unordered_map<int, std::vector<int>> KmerHash;
	typedef std::shared_ptr<SetNode<int>> SetPtr;
	typedef std::unordered_map<int, SetPtr> ClusterHash;
	typedef std::unordered_map<int, std::string> IdToHeader;

	void makeFastaHash(FastaSequences& seqs, FastaHash& hash, IdToHeader& seqsEnum);
	void extractKmers(FastaHash& sequences, KmerHash& kmerHash, int kmerLen);
	void clusterSeqs(KmerHash& kmerHash, FastaHash& fastaHash, ClusterHash& clusters, 
					int kmerLen, int nMissmatch);
	void outputClusters(FastaSet& output);

	int _nMissmatches;
	int _kmerSize;
	
	FastaHash _fastaHash;
	KmerHash _kmerHash;
	ClusterHash _clusterHash;
	IdToHeader _seqEnum;
};

#endif
