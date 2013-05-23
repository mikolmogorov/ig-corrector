all: fasta_clusta

fasta_clusta: fasta_clusta.o fasta.o
	g++ -o fasta_clusta fasta_clusta.o fasta.o

fasta_clusta.o: fasta.h fasta_clusta.cpp disjoint_set.h edit_distance.h
	g++ -c fasta_clusta.cpp -std=c++11 -O2

fasta.o: fasta.h fasta.cpp
	g++ -c fasta.cpp -O2

