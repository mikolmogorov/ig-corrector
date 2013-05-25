#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <list>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include "fasta.h"
#include "disjoint_set.h"

typedef std::vector<std::list<FastaRecord>> FastaSet;

struct pair_hash 
{
    size_t operator() (const std::pair<int,int>& v) const 
	{
        return v.first * 31 + v.second;
    }
};

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
	void splitCliques();
	void spltCluster(std::unordered_set<int> vertex, std::list<std::unordered_set<int>>& out);

	int _nMissmatches;
	int _kmerSize;
	
	FastaHash _fastaHash;
	KmerHash _kmerHash;
	ClusterHash _clusterHash;
	IdToHeader _seqEnum;

	std::unordered_set<std::pair<int, int>, pair_hash> _adjacent;
};



#endif
