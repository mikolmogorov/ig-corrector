#include "clustering.h"
#include <utility>
#include <iostream>
#include <algorithm>
#include <set>
#include "fasta.h"
#include "disjoint_set.h"
#include "edit_distance.h"
#include "bron_kerbosch.h"

void Clusterisator::doJob(FastaSequences& seqs, FastaSet& output, int kmerSize, int nMissmatches)
{
	_fastaHash.clear();
	_kmerHash.clear();
	_clusterHash.clear();
	_seqEnum.clear();
	_kmerSize = kmerSize;
	_nMissmatches = nMissmatches;

	this->makeFastaHash(seqs, _fastaHash, _seqEnum);
	for (auto &itHash : _fastaHash)
	{
		_clusterHash[itHash.first] = SetPtr(makeNode(itHash.first));
	}

	this->extractKmers(_fastaHash, _kmerHash, _kmerSize);
	this->clusterSeqs(_kmerHash, _fastaHash, _clusterHash, _kmerSize, _nMissmatches);
	//this->splitCliques();
	this->outputClusters(output);
}

void Clusterisator::extractKmers(FastaHash& sequences, KmerHash& kmerHash, int kmerLen)
{
	int kmerCounter = 0;
	std::unordered_map<std::string, int> kmerNumbers;

	for (auto &itSeq : sequences)
	{
		for (size_t k = 0; k < itSeq.second.length() - kmerLen; ++k)
		{
			std::string kmer(itSeq.second, k, kmerLen);
			auto itFound = kmerNumbers.find(kmer);
			int kmerNo = 0;
			if (itFound == kmerNumbers.end())
			{
				kmerNumbers.insert(itFound, std::make_pair(kmer, kmerCounter));
				kmerNo = kmerCounter++;
			}
			else
			{
				kmerNo = itFound->second;
			}
			kmerHash[itSeq.first].push_back(kmerNo);
		}
		std::sort(kmerHash[itSeq.first].begin(), kmerHash[itSeq.first].end());
	}
}

void Clusterisator::clusterSeqs(KmerHash& kmerHash, FastaHash& fastaHash, 
								ClusterHash& clusters, int kmerLen, int nMissmatch)
{
	int matchedReads = 0;
	int old_precent = 0;
	int nSeqs = kmerHash.size() * kmerHash.size() / 2;
	int count = 0;

	for (auto firstSeq = kmerHash.begin(); firstSeq != kmerHash.end(); ++firstSeq)
	{
		for (auto secondSeq = firstSeq; secondSeq != kmerHash.end(); ++secondSeq)
		{
			if (firstSeq == secondSeq) continue;

			int matchCount = 0;
			int j = 0, k = 0;
			while (j < firstSeq->second.size() && k < secondSeq->second.size())
			{
				if (firstSeq->second[j] == secondSeq->second[k])
				{
					++matchCount; ++j; ++k;
				}
				else if (firstSeq->second[j] < secondSeq->second[k])
				{
					++j;
				}
				else
				{
					++k;
				}
			}

			int firstKmers = fastaHash[firstSeq->first].length() - kmerLen + 1;
			int secondKmers = fastaHash[secondSeq->first].length() - kmerLen + 1;
			int minSimilarFst = firstKmers - nMissmatch * (2 * kmerLen - 1);
			int minSimilarSnd = secondKmers - nMissmatch * (2 * kmerLen - 1);

			if (matchCount >= std::max(minSimilarFst, minSimilarSnd))
			{
				int eDist = edit_distance(fastaHash[firstSeq->first], fastaHash[secondSeq->first]);
				if (eDist <= nMissmatch)
				{
					++matchedReads;
					auto set1 = findSet(clusters[firstSeq->first].get());
					auto set2 = findSet(clusters[secondSeq->first].get());
					if (set1 != set2)
					{
						unionSet(set1, set2);
					}
					_adjacent.insert(std::make_pair(firstSeq->first, secondSeq->first));
					_adjacent.insert(std::make_pair(secondSeq->first, firstSeq->first));
				}
			}

			++count;
			int percent = float(count) / nSeqs * 100;
			if (percent > old_precent)
			{
				std::cerr << percent << " ";
				old_precent = percent;
				if (percent % 10 == 0) std::cerr << std::endl;
			}
		}
	}
	std::cerr << std::endl;
}

void Clusterisator::makeFastaHash(FastaSequences& seqs, FastaHash& hash, IdToHeader& seqsEnum)
{
	long long counter = 0;
	for (auto &record : seqs)
	{
		hash[counter] = record.sequence;
		seqsEnum[counter] = record.description;
		++counter;
	}
}

void Clusterisator::outputClusters(FastaSet& output)
{
	std::unordered_map<int, std::unordered_set<int>> preClusters;
	std::list<std::unordered_set<int>> clusters;

	for (auto &item : _clusterHash)
	{
		auto parent = findSet(item.second.get());
		preClusters[parent->data].insert(item.first);
	}

	for (auto &cl : preClusters)
	{
		this->spltCluster(cl.second, clusters);
	}

	for (auto &itClust : clusters)
	{
		output.push_back(std::list<FastaRecord>());
		for (int seqId : itClust)
		{
			output.back().push_back(FastaRecord(_fastaHash[seqId], _seqEnum[seqId], seqId));
		}
	}
}

void Clusterisator::spltCluster(std::unordered_set<int> vertex, std::list<std::unordered_set<int>>& out)
{
	std::cerr << "== cluster of " << vertex.size() << std::endl;
	//while (true)
	while (!vertex.empty())
	{
		/*
		if (vertex.size() < 3)
		{
			if (!vertex.empty()) out.push_back(vertex);
			return;
		}*/

		Graph g;

		for (auto v1 = vertex.begin(); v1 != vertex.end(); ++v1)
		{
			for (auto v2 = v1; v2 != vertex.end(); ++v2)
			{
				if (v1 == v2) continue;
				if (_adjacent.find(std::make_pair(*v1, *v2)) != _adjacent.end())
				{
					g.add_edge(*v1, *v2);
					g.add_edge(*v2, *v1);
				}
			}
		}
		std::vector<Graph::vertex_set> vvs;
		g.find_cliques(vvs);

		if (vvs.empty())
		{
			std::cerr << "\tno more cliques, " << vertex.size() << " left" << std::endl;
			//out.push_back(vertex);
			for (auto elem : vertex)
			{
				out.push_back(std::unordered_set<int>());
				out.back().insert(elem);
			}
			return;
		}

		Graph::vertex_set* maxClique = &vvs.front();
		for (Graph::vertex_set &vs : vvs)
		{
			if (vs.size() > maxClique->size())
			{
				maxClique = &vs;
			}
		}
		
		std::cerr << "\tclique of " << maxClique->size() << std::endl;
		out.push_back(std::unordered_set<int>(maxClique->begin(), maxClique->end()));

		for (auto elem : *maxClique)
		{
			vertex.erase(elem);
		}
	}
}

/*
void Clusterisator::splitCliques()
{
	std::unordered_map<int, std::unordered_set<int>> clusters;
	for (auto &item : _clusterHash)
	{
		auto parent = findSet(item.second.get());
		clusters[parent->data].insert(item.first);
	}

	for (auto &itClust : clusters)
	{
		spltCluster(itClust.second);
	}
}*/
