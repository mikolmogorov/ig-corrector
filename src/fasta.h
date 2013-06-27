#ifndef _FASTA_READER_H_
#define _FASTA_READER_H_

#include <string>
#include <vector>
#include <fstream>
#include <iterator>

struct FastaRecord
{
	FastaRecord() {}
	FastaRecord(const std::string & sequence, const std::string & description, size_t id):
		sequence(sequence), description(description), seqId(id)
	{
	}
	
	std::string sequence;
	std::string description;		
	size_t seqId;
};

typedef std::vector<FastaRecord> FastaSequences;
typedef std::vector<FastaSequences> FastaSet;

class FastaReader
{
public:
	explicit FastaReader(std::istream& stream): 
		inputStream_(stream)
	{}
	size_t 	GetSequences(std::vector<FastaRecord>& record);
	bool 	IsOk() const;
private:
	void ValidateSequence(std::string & sequence);
	void ValidateHeader(std::string & header);

	std::istream& inputStream_;
};

class FastaWriter
{
public:
	static void WriteSequence(std::ostream& out, 
								const std::string & header, 
								const std::string & sequence)
	{
		out << ">" << header << std::endl;
		for(size_t i = 0; i < sequence.size(); i += 80)
		{
			size_t j = std::min(i + 80, sequence.size());
			std::copy(sequence.begin() + i, sequence.begin() + j, std::ostream_iterator<char>(out));
			out << std::endl;
		}
	}
};

inline void writeFastaSet(const FastaSet& set, std::ostream& out)
{
	int counter = 0;
	for (FastaSet::const_iterator itCluster = set.begin(); itCluster != set.end(); ++itCluster)
	{
		out << "=Cluster_" << counter << "_" << itCluster->size() << std::endl;
		for (FastaSequences::const_iterator itSeq = itCluster->begin(); itSeq != itCluster->end(); ++itSeq)
		{
			FastaWriter::WriteSequence(out, itSeq->description, itSeq->sequence);
		}
		++counter;
	}
}

#endif 
