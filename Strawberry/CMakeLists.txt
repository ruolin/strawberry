# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Strawberry )


cmake_minimum_required(VERSION 2.8.12)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
   add_compile_options(-Wall -std=c++11)
elseif(COMPILER_SUPPORTS_CXX0X)
   add_complie_options(-Wall -DDEBUG -std=c++0x)
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package


# Creating entries for target: Strawberry
# ############################

#
set(LEMON_INCLUDE lemon/include)
set(BAM_INCLUDE /usr/local/include/bam)
link_directories(/usr/local/lib)
include_directories(include )
add_library(strawlib common.cpp fasta.cpp contig.cpp read.cpp gff.cpp logger.cpp qp.cpp)

add_library(alignments.o OBJECT alignments.cpp include/alignments.h)
target_include_directories(alignments.o PRIVATE ${LEMON_INCLUDE})
add_library(assembly.o OBJECT assembly.cpp include/assembly.h include/common.h)
target_include_directories(assembly.o PRIVATE ${LEMON_INCLUDE})

add_executable(Strawberry Strawberry.cpp  $<TARGET_OBJECTS:alignments.o> $<TARGET_OBJECTS:assembly.o> ) 
#add_to_cached_list( CGAL_EXECUTABLE_TARGETS Strawberry )
# Link the executable to CGAL and third-party libraries
target_link_libraries(Strawberry strawlib bam z pthread ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

