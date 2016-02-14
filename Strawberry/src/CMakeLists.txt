link_directories(
  /usr/local/lib
  ${BAM_LIBRARY}
)

include_directories(
  ${PROJECT_SOURCE_DIR}/include
 # ${PROJECT_SOURCE_DIR}/include/eigen
 # ${PROJECT_SOURCE_DIR}/include/lemon
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
)

#add_library(strawlib common.cpp fasta.cpp contig.cpp read.cpp gff.cpp logger.cpp)
#add_library(estimate.o OBJECT estimate.cpp assembly.h)
#target_include_directories(estimate.o PRIVATE)
#add_library(alignments.o OBJECT alignments.cpp alignments.h)
#target_include_directories(alignments.o PRIVATE)
#add_library(assembly.o OBJECT assembly.cpp assembly.h common.h)
#target_include_directories(assembly.o PRIVATE) 
add_executable(strawberry ${STRAWBERRY_MAIN_SRCS})
target_link_libraries(strawberry bam z pthread )

