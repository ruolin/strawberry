CXX = g++
CXXFLAGS =	-O3 -g -std=c++11 -Wall 
INCLUDES = -I /usr/local/include
LIBS = -L /usr/local/lib/ -l bam -l z -l pthread
LDLIBS= -lbam
MAIN_TARGET =	Strawberry

read.o: read.cpp read.h
genome.o : genome.cpp genome.h
	${CXX} ${CXXFLAGS} -c genome.cpp
GBase.o : GBase.cpp GBase.h
	${CXX} ${CXXFLAGS} -c GBase.cpp

Strawberry.o : Strawberry.cpp GBase.h genome.h read.h
	${CXX} ${CXXFLAGS} -c Strawberry.cpp
	
$(MAIN_TARGET):	Strawberry.o GBase.o genome.o read.o
	$(CXX) ${CXXFLAGS} Strawberry.o GBase.o genome.o read.o -o $(MAIN_TARGET) $(INCLUDES)  $(LIBS) 

all:	$(MAIN_TARGET)

clean:
	rm -f $(MAIN) $(MAIN_TARGET)
	rm -f *.o
