all: read_cluster

read_cluster: read_cluster.o fasta.o
	g++ -o read_cluster read_cluster.o fasta.o

read_cluster.o: fasta.h read_cluster.cpp disjoint_set.h edit_distance.h
	g++ -c read_cluster.cpp -std=c++11 -O2

fasta.o: fasta.h fasta.cpp
	g++ -c fasta.cpp -O2

