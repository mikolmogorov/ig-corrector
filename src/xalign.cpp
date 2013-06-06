#include <iostream>
#include "xalign_impl.h"
#include <getopt.h>
#include "fasta.h"

bool parseArgs(int argc, char** argv, int& trhld, int& match, int& missmatch, 
				int& ins, int& del, std::string& querry, std::string& filename);

int main(int argc, char* argv[])
{
	int threshold;
	int matchScore = 2;
	int missmatchScore = -2;
	int delScore = -1;
	int insScore = -1;
	std::string querry;
	std::string fileName;

	if (!parseArgs(argc, argv, threshold, matchScore, missmatchScore, insScore, delScore, querry, fileName))
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
	//std::cout << matchScore << " " << missmatchScore;

	Xaligner aligner(matchScore, missmatchScore, insScore, delScore);
	aligner.setQuerry(querry);
	for (auto seq : seqs)
	{
		std::vector<AlignInfo> result;
		aligner.align(seq.sequence, threshold, result);
		std::cout << ">" << seq.description << " ";
		for (AlignInfo& aln : result)
		{
			std::cout << "(" << aln.start << " " << aln.end << " " << aln.score << ") ";
		}
		std::cout << std::endl;
	}
/*	
	Xaligner aligner;
	aligner.setQuerry("YYC");
	std::vector<AlignInfo> result;
	aligner.align("TATTACTGTACAAATCCCACAGTCCGGGACTACTTGGGCAAAGGGACCACGGTCACCGTCTCCTCAGGCTCGAGTGCGTCTACAA", 
					15, result);
	for (auto aln : result)
	{
		std::cout << aln.start << " " << aln.end << " " << aln.score << std::endl;
	}*/
	return 0;
}

bool parseArgs(int argc, char** argv, int& trhld, int& match, int& missmatch, 
				int& ins, int& del, std::string& querry, std::string& filename)
{
	auto printUsage = []()
	{
		std::cerr 	<< "\nUSAGE: xalign [-m match_score] [-x missmatch_score] "
					<< "[-i insert_score]\n\t\t[-d delete_score] -t threshold -q querry reads_file\n"
					<< "Default values: m = 2; x = -2; d = -1; i = -1\n"
					<< "If reads_file is not set, reading from standard input\n\n";
	};
	const char* optString = "m:x:i:d:t:q:?";
	int opt = 0;
	bool thresholdSet = false;
	while( (opt = getopt(argc, argv, optString)) != -1 )
	{
		switch(opt)
		{
		case 'm':
			match = atoi(optarg);
			break;
		case 'x':
			missmatch = atoi(optarg);
			break;
		case 'i':
			ins = atoi(optarg);
			break;
		case 'd':
			del = atoi(optarg);
			break;
		case 't':
			trhld = atoi(optarg);
			thresholdSet = true;
			break;
		case 'q':
			querry = optarg;
			break;
		case '?':
			printUsage();
			return false;
			break;
		}
	}
	if (!thresholdSet || querry.empty())
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
