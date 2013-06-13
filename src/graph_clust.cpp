#include <iostream>
#include <stdexcept>
#include <getopt.h>
#include "fasta.h"
#include "graph_impl.h"

bool parseArgs(int argc, char** argv, int& kmerSize, int& nMissmatch, bool& quiet, std::string& filename);

int main(int argc, char* argv[])
{
	int kmer = 0;
	int threshold = 0;
	bool quiet = false;
	std::string fileName;

	if (!parseArgs(argc, argv, kmer, threshold, quiet, fileName))
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
			std::cerr << "error opening " << fileName << std::endl;
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
		std::cerr << e.what() << std::endl;
		return 1;
	}

	FastaSet outClusters;
	Clusterisator clusterisator;
	clusterisator.doJob(seqs, outClusters, kmer, threshold, !quiet);
	writeFastaSet(outClusters, std::cout);
	return 0;
}

bool parseArgs(int argc, char** argv, int& kmerSize, int& nMissmatch, bool& quiet, std::string& filename)
{
	auto printUsage = []()
	{
		std::cerr 	<< "\nUSAGE: graph_clust -k kmer_size -m num_missmatch [-q] reads_file\n"
					<< "If reads_file is not set, reading from standard input\n\n";
	};
	const char* optString = "k:m:q?";
	int opt = 0;
	bool kmerSet = false;
	bool missSet = false;
	quiet = false;
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
		case 'q':
			quiet = true;
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
