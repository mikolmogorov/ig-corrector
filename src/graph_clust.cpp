#include <iostream>
#include <stdexcept>
#include <getopt.h>
#include "fasta.h"
#include "graph_impl.h"

//void outputClusters(FastaSet& clusters);

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

	FastaSet outClusters;
	Clusterisator clusterisator;
	clusterisator.doJob(seqs, outClusters, KMER, THRESHOLD);
	writeFastaSet(outClusters, std::cout);
	//outputClusters(outClusters);
	return 0;
}

bool parseArgs(int argc, char** argv, int& kmerSize, int& nMissmatch, std::string& filename)
{
	auto printUsage = []()
	{
		std::cerr 	<< "\nUSAGE: graph_clust -k kmer_size -m num_missmatch reads_file\n"
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

/*
void outputClusters(FastaSet& clusters)
{
	int counter = 0;
	for (auto &clust : clusters)
	{
		std::cout << "Cluster" << counter << "_" << clust.size() << std::endl;
		for (auto &itSeq : clust)
		{
			std::cout << ">" << itSeq.description << std::endl << itSeq.sequence << std::endl;
		}
		++counter;
		std::cout << std::endl;
	}
}*/
