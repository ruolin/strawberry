CXX = g++
CXXFLAGS =	-O3 -g -std=c++11 -Wall -DDEBUG
INCLUDES = -I /usr/local/include
LIBS = -L /usr/local/lib/ -l bam -l z -l pthread
LDLIBS= -lbam
MAIN_TARGET =	Strawberry

fasta.o: fasta.h common.h fasta.cpp
	${CXX} ${CXXFLAGS} -c fasta.cpp
contig.o: contig.h contig.cpp
	${CXX} ${CXXFLAGS} -c contig.cpp
read.o: read.h contig.h read.cpp 
	${CXX} ${CXXFLAGS} -c read.cpp
gff.o : gff.h gff.cpp 
	${CXX} ${CXXFLAGS} -c gff.cpp
GBase.o : GBase.h GBase.cpp 
	${CXX} ${CXXFLAGS} -c GBase.cpp

Strawberry.o : Strawberry.cpp GBase.h gff.h read.h fasta.h
	${CXX} ${CXXFLAGS} -c Strawberry.cpp
	
$(MAIN_TARGET):	Strawberry.o GBase.o gff.o read.o contig.o fasta.o
	$(CXX) ${CXXFLAGS} Strawberry.o GBase.o gff.o read.o contig.o fasta.o -o $(MAIN_TARGET) $(INCLUDES)  $(LIBS) 

all:	$(MAIN_TARGET)

clean:
	rm -f $(MAIN) $(MAIN_TARGET)
	rm -f *.o
