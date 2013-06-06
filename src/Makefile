all: graph_clust hierarchial_clust xalign

graph_clust: graph_clust.o fasta.o graph_impl.o
	g++ -o graph_clust graph_clust.o fasta.o graph_impl.o 

graph_clust.o: fasta.h graph_clust.cpp graph_impl.h
	g++ -c graph_clust.cpp -std=c++11 -O2

fasta.o: fasta.h fasta.cpp
	g++ -c fasta.cpp -O2

graph_impl.o: graph_impl.cpp graph_impl.h fasta.h common.h
	g++ -c graph_impl.cpp -std=c++11 -O2

hierarchial_clust: hierarchial_clust.o fasta.o hierarchial_impl.o
	g++ -o hierarchial_clust hierarchial_clust.o fasta.o hierarchial_impl.o

hierarchial_clust.o: fasta.h hierarchial_impl.h hierarchial_clust.cpp
	g++ -c hierarchial_clust.cpp -std=c++11 -O2

hierarchial_impl.o: fasta.h edit_distance.h hierarchial_impl.cpp hierarchial_impl.h common.h
	g++ -c hierarchial_impl.cpp -std=c++11 -O2

xalign: xalign.o xalign_impl.o fasta.o
	g++ -o xalign xalign.o xalign_impl.o fasta.o

xalign.o: xalign_impl.h xalign.cpp fasta.h
	g++ -c xalign.cpp -std=c++11 -O2

xalign_impl.o: xalign_impl.cpp xalign_impl.h
	g++ -c xalign_impl.cpp -std=c++11 -O2
