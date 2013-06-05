#include "fasta.h"
#include "hierarchial_impl.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[])
{
	std::ifstream fin(argv[1]);
	FastaReader reader(fin);
	FastaSequences seqs;
	FastaSet clusters;
	try
	{
		reader.GetSequences(seqs);
	}
	catch (std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return 1;
	}
	HierarchialClust clusterizator;
	clusterizator.cluster(seqs, clusters, 4.0f);
	writeFastaSet(clusters, std::cout);

	return 0;
}
