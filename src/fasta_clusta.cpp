#include <unordered_map>
#include <vector>
#include <string>
#include <list>
#include <utility>
#include <iostream>
#include <algorithm>
#include <memory>
#include <set>
#include <getopt.h>
#include "fasta.h"
#include "disjoint_set.h"
#include "edit_distance.h"
//#include <cassert>

typedef std::vector<FastaRecord> FastaSequences;
typedef std::unordered_map<std::string, std::string> FastaHash;
typedef std::unordered_map<std::string, std::vector<int>> KmerHash;
typedef std::shared_ptr<SetNode<std::string>> SetPtr;
typedef std::unordered_map<std::string, SetPtr> ClusterHash;

void extractKmers(FastaSequences& sequences, KmerHash& hash, int kmerLen);
void clusterSeqs(KmerHash& kmerHash, FastaHash& fastaHash, ClusterHash& clusters, int kmerLen, int nMissmatch);
void makeFastaHash(FastaSequences& seqs, FastaHash& hash);
void outputClusters(FastaHash& fastaHash, ClusterHash& clusterHash);

bool parseArgs(int argc, char** argv, int& kmerSize, int& nMissmatch, std::string& filename);

int main(int argc, char* argv[])
{
	int KMER = 0;
	int THRESHOLD = 0;
	std::string fileName;

	if (!parseArgs(argc, argv, KMER, THRESHOLD, fileName))
	{
		return 1;
	}

	FastaSequences seqs;
	FastaHash seqHash;
	KmerHash kmerHash;
	ClusterHash clusterHash;
	
	std::istream* stream = nullptr;
	std::ifstream fstream;
	if (fileName.empty())
	{
		stream = &std::cin;
	}
	else
	{
		fstream.open(fileName);
		if (!fstream.is_open())
		{
			std::cout << "error opening " << fileName << std::endl;
			return 1;
		}
		stream = &fstream;
	}

	FastaReader reader(*stream);
	try
	{
		reader.GetSequences(seqs);
	}
	catch (std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return 1;
	}

	for (auto &seq : seqs)
	{
		clusterHash[seq.description] = SetPtr(makeNode(seq.description));
	}

	makeFastaHash(seqs, seqHash);
	extractKmers(seqs, kmerHash, KMER);
	clusterSeqs(kmerHash, seqHash, clusterHash, KMER, THRESHOLD);
	outputClusters(seqHash, clusterHash);
	return 0;
}

bool parseArgs(int argc, char** argv, int& kmerSize, int& nMissmatch, std::string& filename)
{
	auto printUsage = []()
	{
		std::cerr 	<< "\nUSAGE: fasta_clusta -k kmer_size -m num_missmatch reads_file\n"
					<< "If reads_file is not set, reading from standard input\n\n";
	};
	const char* optString = "k:m:?";
	int opt = 0;
	bool kmerSet = false;
	bool missSet = false;
	while( (opt = getopt(argc, argv, optString)) != -1 )
	{
		switch(opt)
		{
		case 'k':
			kmerSize = atoi(optarg);
			kmerSet = true;
			break;
		case 'm':
			nMissmatch = atoi(optarg);
			missSet = true;
			break;
		case '?':
			printUsage();
			return false;
			break;
		}
	}
	if (!kmerSet || !missSet)
	{
		printUsage();
		return false;
	}
	if (argc - optind > 0)
	{
		filename = *(argv + optind);
	}
	return true;
}

void extractKmers(FastaSequences& sequences, KmerHash& hash, int kmerLen)
{
	int kmerCounter = 0;
	std::unordered_map<std::string, int> kmerNumbers;

	for (FastaRecord& record : sequences)
	{
		//assert(hash.find(record.description) == hash.end());
		for (size_t k = 0; k < record.sequence.length() - kmerLen; ++k)
		{
			std::string kmer(record.sequence, k, kmerLen);
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
			hash[record.description].push_back(kmerNo);
		}
		std::sort(hash[record.description].begin(), hash[record.description].end());
	}
}

void clusterSeqs(KmerHash& kmerHash, FastaHash& fastaHash, ClusterHash& clusters, int kmerLen, int nMissmatch)
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

void makeFastaHash(FastaSequences& seqs, FastaHash& hash)
{
	for (auto &record : seqs)
	{
		hash[record.description] = record.sequence;
	}
}

void outputClusters(FastaHash& fastaHash, ClusterHash& clusterHash)
{
	std::unordered_map<std::string, std::list<std::string>> clusters;
	for (auto &item : clusterHash)
	{
		auto parent = findSet(item.second.get());
		clusters[parent->data].push_back(item.first);
	}
	int counter = 0;
	for (auto &clust : clusters)
	{
		std::cout << "Cluster" << counter << "_" << clust.second.size() << std::endl;
		for (auto &item : clust.second)
		{
			std::cout << ">" << item << std::endl << fastaHash[item] << std::endl;
		}
		++counter;
		std::cout << std::endl;
	}
}
