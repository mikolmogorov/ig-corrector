CC=g++		#change with your compiler name
CFLAGS=-c -Wall -std=c++0x -O2
LDFLAGS=

#change on this, if you are running troubles with clang
#CFLAGS=-c -Wall -std=c++0x -stdlib=libc++ -O2 	
#LDFLAGS=-lc++

all: graph_clust hierarchial_clust xalign

graph_clust: graph_clust.o fasta.o graph_impl.o
	$(CC) $(LDFLAGS) -o graph_clust graph_clust.o fasta.o graph_impl.o

graph_clust.o: fasta.h graph_clust.cpp graph_impl.h
	$(CC) $(CFLAGS) graph_clust.cpp

fasta.o: fasta.h fasta.cpp
	$(CC) $(CFLAGS) fasta.cpp

graph_impl.o: graph_impl.cpp graph_impl.h fasta.h common.h
	$(CC) $(CFLAGS) graph_impl.cpp

hierarchial_clust: hierarchial_clust.o fasta.o hierarchial_impl.o
	$(CC) $(LDFLAGS) -o hierarchial_clust hierarchial_clust.o fasta.o hierarchial_impl.o

hierarchial_clust.o: fasta.h hierarchial_impl.h hierarchial_clust.cpp
	$(CC) $(CFLAGS) hierarchial_clust.cpp

hierarchial_impl.o: fasta.h edit_distance.h hierarchial_impl.cpp hierarchial_impl.h common.h
	$(CC) $(CFLAGS) hierarchial_impl.cpp

xalign: xalign.o xalign_impl.o fasta.o
	$(CC) $(LDFLAGS) -o xalign xalign.o xalign_impl.o fasta.o

xalign.o: xalign_impl.h xalign.cpp fasta.h
	$(CC) $(CFLAGS) xalign.cpp

xalign_impl.o: xalign_impl.cpp xalign_impl.h
	$(CC) $(CFLAGS) xalign_impl.cpp

clean:
	rm -rf *.o graph_clust hierarchial_clust xalign
