#include "graph_impl.h"
#include <utility>
#include <iostream>
#include <algorithm>
#include <set>
#include "fasta.h"
#include "disjoint_set.h"
#include "edit_distance.h"

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
		clusters.push_back(cl.second);
	}

	for (auto &itClust : clusters)
	{
		output.push_back(std::vector<FastaRecord>());
		for (int seqId : itClust)
		{
			output.back().push_back(FastaRecord(_fastaHash[seqId], _seqEnum[seqId], seqId));
		}
	}
}
