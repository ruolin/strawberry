CXX = g++
CXXFLAGS =	-O3 -g -std=c++11 -Wall -DDEBUG
INCLUDES = -I/usr/local/include/bam 
LEMON_INCLUDE = -I lemon/include
#LEMON_LIB = -L lemon/lib
LIBS = -L/usr/local/lib/ -l bam -l z -l pthread 
LDLIBS= -lbam
MAIN_TARGET = Strawberry


common.o: common.h
	${CXX} ${CXXFLAGS} -c common.cpp
fasta.o: fasta.h common.h fasta.cpp
	${CXX} ${CXXFLAGS} -c fasta.cpp
contig.o: contig.h contig.cpp common.h
	${CXX} ${CXXFLAGS} -c contig.cpp
read.o: read.h contig.h read.cpp 
	${CXX} ${CXXFLAGS} -c read.cpp
gff.o : gff.h gff.cpp contig.h
	${CXX} ${CXXFLAGS} -c gff.cpp
alignments.o: alignments.h alignments.cpp 
	${CXX} ${CXXFLAGS} -c alignments.cpp $(LEMON_INCLUDE)
assembly.o: assembly.h assembly.cpp
	${CXX} ${CXXFLAGS} -c assembly.cpp $(LEMON_INCLUDE)
logger.o: logger.cpp logger.hpp
	${CXX} ${CXXFLAGS} -c logger.cpp

Strawberry.o : Strawberry.cpp gff.h read.h fasta.h alignments.h logger.hpp common.h assembly.o
	${CXX} ${CXXFLAGS} -c Strawberry.cpp 
	
$(MAIN_TARGET):	Strawberry.o gff.o read.o contig.o fasta.o common.o alignments.o logger.o assembly.o
	$(CXX) ${CXXFLAGS} Strawberry.o gff.o read.o contig.o fasta.o common.o alignments.o logger.o assembly.o -o $(MAIN_TARGET) $(INCLUDES) $(LIBS) 

all:	$(MAIN_TARGET)

clean:
	rm -f $(MAIN) $(MAIN_TARGET)
	rm -f *.o
