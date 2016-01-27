cmake_minimum_required(VERSION 2.8.12)

project( Strawberry )

set(CPACK_PACKAGE_VERSION "0.7.2")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "7")
set(CPACK_PACKAGE_VERSION_PATCH "2")
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "Iowa State University")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Strawberry - Fast and accurate isoform reconstraction and quantification using RNA-seq")
set(CPACK_PACKAGE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-Source")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
   add_compile_options(-Wall -g -O3 -DDEBUG -std=c++11)
elseif(COMPILER_SUPPORTS_CXX0X)
   add_complie_options(-Wall -g -DDEBUG -std=c++0x)
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()


# include for local directory

# include for local package


# Creating entries for target: Strawberry
# ############################

#

set(BAM_INCLUDE /usr/local/include/bam)
link_directories(/usr/local/lib)
include_directories(include )
add_library(strawlib common.cpp fasta.cpp contig.cpp read.cpp gff.cpp logger.cpp)
add_library(estimate.o OBJECT estimate.cpp include/assembly.h)
target_include_directories(estimate.o PRIVATE)
add_library(alignments.o OBJECT alignments.cpp include/alignments.h)
target_include_directories(alignments.o PRIVATE)
add_library(assembly.o OBJECT assembly.cpp include/assembly.h include/common.h)
target_include_directories(assembly.o PRIVATE)

add_executable(Strawberry Strawberry.cpp  $<TARGET_OBJECTS:alignments.o> $<TARGET_OBJECTS:assembly.o> $<TARGET_OBJECTS:estimate.o>) 
target_link_libraries(Strawberry strawlib bam z pthread )

