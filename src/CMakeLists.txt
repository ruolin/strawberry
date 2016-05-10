#link_directories(
  #/usr/local/lib
  #${BAM_LIBRARY}
#)

add_library(bam STATIC IMPORTED)
set_property(TARGET bam PROPERTY IMPORTED_LOCATION ${BAM_LIBRARY})

include_directories(

  #${PROJECT_SOURCE_DIR}/include/eigen
  #${PROJECT_SOURCE_DIR}/include/sam
  ${PROJECT_SOURCE_DIR}/include
)

set( STRAWBERRY_MAIN_SRCS
common.cpp
fasta.cpp
contig.cpp
read.cpp
gff.cpp
logger.cpp
estimate.cpp
alignments.cpp
assembly.cpp
Strawberry.cpp
kmer.cpp
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
#add_library(strawlib common.cpp fasta.cpp contig.cpp read.cpp gff.cpp logger.cpp)
#add_library(estimate.o OBJECT estimate.cpp assembly.h)
#target_include_directories(estimate.o PRIVATE)
#add_library(alignments.o OBJECT alignments.cpp alignments.h)
#target_include_directories(alignments.o PRIVATE)
#add_library(assembly.o OBJECT assembly.cpp assembly.h common.h)
#target_include_directories(assembly.o PRIVATE) 

add_executable(strawberry ${STRAWBERRY_MAIN_SRCS})
#target_link_libraries(bam z)
target_link_libraries(strawberry bam z pthread )

install(TARGETS strawberry 
			    RUNTIME DESTINATION bin
)
