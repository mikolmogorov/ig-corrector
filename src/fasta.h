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
		description(description), seqId(id), sequence(sequence)
	{
	}
	
	std::string sequence;
	std::string description;		
	size_t seqId;
};

class FastaReader
{
public:
	explicit FastaReader(const std::string & fileName): 
		inputStream_(fileName.c_str()),
		fileName_(fileName) {}

	size_t 	GetSequences(std::vector<FastaRecord>& record);
	bool 	IsOk() const;
private:
	void ValidateSequence(std::string & sequence);
	void ValidateHeader(std::string & header);

	std::ifstream inputStream_;
	std::string fileName_;
};

class FastaWriter
{
public:
	static void WriteSequence(const std::string & fileName, 
								const std::string & header, 
								const std::string & sequence)
	{
		std::ofstream out(fileName.c_str());
		out << ">" << header << std::endl;
		for(size_t i = 0; i < sequence.size(); i += 80)
		{
			size_t j = std::min(i + 80, sequence.size());
			std::copy(sequence.begin() + i, sequence.begin() + j, std::ostream_iterator<char>(out));
			out << std::endl;
		}
	}
};

#endif 
