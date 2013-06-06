#include "fasta.h"
#include "hierarchial_impl.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <getopt.h>

bool parseArgs(int argc, char** argv, float& cutoff, bool& quiet, std::string& filename);

int main(int argc, char* argv[])
{
	float cutoff;
	bool quiet;
	std::string fileName;
	if (!parseArgs(argc, argv, cutoff, quiet, fileName))
	{
		return 1;
	}

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

	FastaSequences seqs;
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

	FastaSet clusters;
	HierarchialClust clusterizator;
	clusterizator.cluster(seqs, clusters, 4.0f);
	writeFastaSet(clusters, std::cout);

	return 0;
}

bool parseArgs(int argc, char** argv, float& cutoff, bool& quiet, std::string& filename)
{
	auto printUsage = []()
	{
		std::cerr 	<< "\nUSAGE: hierarchial_clust -c cutoff [-q] reads_file\n"
					<< "If reads_file is not set, reading from standard input\n\n";
	};
	const char* optString = "c:q?";
	int opt = 0;
	bool cutoffSet = false;
	quiet = false;
	while( (opt = getopt(argc, argv, optString)) != -1 )
	{
		switch(opt)
		{
		case 'c':
			cutoff = atof(optarg);
		 	cutoffSet = true;
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
	if (!cutoffSet)
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
