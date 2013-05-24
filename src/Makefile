all: fasta_clusta

fasta_clusta: fasta_clusta.o fasta.o clustering.o
	g++ -o fasta_clusta fasta_clusta.o fasta.o clustering.o

fasta_clusta.o: fasta.h fasta_clusta.cpp clustering.h
	g++ -c fasta_clusta.cpp -std=c++11 -O2

fasta.o: fasta.h fasta.cpp
	g++ -c fasta.cpp -O2

clustering.o: clustering.cpp clustering.h fasta.h
	g++ -c clustering.cpp -std=c++11 -O2
