#link_directories(
  #/usr/local/lib
  #${BAM_LIBRARY}
#)

add_library(bam STATIC IMPORTED)
set_property(TARGET bam PROPERTY IMPORTED_LOCATION ${BAM_LIBRARY})

include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/external/eigen3
  ${PROJECT_SOURCE_DIR}/external/easyloggingpp/src
)

set( STRAWBERRY_MAIN_SRCS
common.cpp
fasta.cpp
contig.cpp
read.cpp
gff.cpp
estimate.cpp
alignments.cpp
assembly.cpp
Strawberry.cpp
${PROJECT_SOURCE_DIR}/external/easyloggingpp/src/easylogging++.cc
bias.cpp 
../include/isoform.h 
../include/utils.h 
../include/kmer.h)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_executable(strawberry ${STRAWBERRY_MAIN_SRCS})
#target_link_libraries(bam z)
target_link_libraries(strawberry bam z pthread )

install(TARGETS strawberry 
			    RUNTIME DESTINATION bin
)
