all: anchors

anchors: anchors.o fasta.o
	g++ -o anchors anchors.o fasta.o

anchors.o: fasta.h anchors.cpp disjoint_set.h edit_distance.h
	g++ -c anchors.cpp -std=c++11 -O2

fasta.o: fasta.h fasta.cpp
	g++ -c fasta.cpp -O2

