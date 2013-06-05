#ifndef HIERARCH_CLUST_H
#define HIERARCH_CLUST_H

#include <list>
#include <map>
#include <unordered_map>
#include "fasta.h"
#include "common.h"

class HierarchialClust
{
public:
	void cluster(const FastaSequences& seqs, FastaSet& output, float cutoff);
	
private:
	typedef std::pair<int, int> KeyPair;
	typedef std::unordered_map<int, std::list<std::string>> Clusters;
	typedef std::unordered_map<std::string, std::string> FastaHash;
	typedef std::unordered_map<KeyPair, float, pair_hash<int>> DistHash;

	bool 	step();
	int 	sequenceDistance(const std::string& seq1, const std::string& seq2);
	float	clusterDist(int clust1, int clust2);
	void	outputClusters(FastaSet& out);
	
	
	DistHash _distances; 
	Clusters _clusters;
	FastaHash _fastaHash; 
	float _cutoff;
};


#endif
