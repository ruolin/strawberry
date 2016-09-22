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
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
#MESSAGE( STATUS "CMAKE_MODULE_DIR:         " ${CMAKE_MODULE_PATH} )
set(CMAKE_CXX_FLAGS "-g -pthread -Ofast -funroll-loops -fPIC -fomit-frame-pointer -Wall -Wno-reorder -Wno-unused-variable -Wreturn-type -Werror=return-type")
if (APPLE)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
else()
    # We're in sane linux world
endif()

message("Using g++ complier: " ${CMAKE_CXX_COMPILER})
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    # If we're on OSX
    if (APPLE AND NOT (GCC_VERSION VERSION_GREATER 4.8.2 OR GCC_VERSION VERSION_EQUAL 4.8.2))
        message(FATAL_ERROR "When building under OSX, ${PROJECT_NAME} requires "
                            "g++ >= 4.8.2")
    elseif (NOT (GCC_VERSION VERSION_GREATER 4.8.2 OR GCC_VERSION VERSION_EQUAL 4.8.2))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.8.2 or greater.")
    endif ()

    set (GCC TRUE)

    # Put complete static linking on hold for the time-being
    # If we're not on OSX, make an attempt to compile everything statically
    set (PTHREAD_LIB "pthread")

    # If we're on Linux (i.e. not OSX) and we're using
    # gcc, then set the -static-libstdc++ flag
    if (NOT APPLE)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")
    endif()

    set (WARNING_IGNORE_FLAGS "${WARNING_IGNORE_FLAGS} -Wno-unused-local-typedefs")
else()
	 message(FATAL_ERROR "Please install g++.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_IGNORE_FLAGS}")

##
# Record this top-level path
##
set (GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Have CMake tell us what it's doing
set (CMAKE_VERBOSE_MAKEFILE true)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)


if(COMPILER_SUPPORTS_CXX11)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()



if (CMAKE_BUILD_TYPE MATCHES Debug)
    message ("Making Debug build")
    set (CMAKE_CXX_FLAGS_DEBUG "-g ${CMAKE_CXX_FLAGS}")
elseif (CMAKE_BUILD_TYPE MATCHES Release)
    message ("Making Release build")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS}")
else ()
    message ("Making Default build type")
endif ()
#message("C++ complier options: ${CMAKE_CXX_FLAGS_RELEASE}")

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

find_package (ZLIB)
if (NOT ZLIB_FOUND)
	message (FATAL_ERROR "zlib must be installed before configuration & building can proceed")
endif()

# include for local directory

# include for local package


# Creating entries for target: Strawberry
# ############################

include(ExternalProject)


#set(BAM_INCLUDE /usr/local/include/bam)


find_library(BAM_LIBRARY 
    NAMES libbam bam bamlib
	HINTS
		${SAMTOOLS_ROOT} ENV SAMTOOLS_ROOT
	PATHS
		#/usr/local/lib
        ${CMAKE_CURRENT_SOURCE_DIR}/exteral/samtools-0.1.19
)


if( BAM_LIBRARY)
	message("-- Found samtools: ${BAM_LIBRARY}")
else()
	message("-- samtools library not found! ")
    message("-- will attempt to fetch and build samtools on version 0.1.19")
ExternalProject_Add(samtools
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
  URL https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools-0.1.19
# INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
# PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools-0.1.19
  BUILD_IN_SOURCE TRUE
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
  INSTALL_COMMAND ""
)

set(BAM_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools-0.1.19/libbam.a)
endif()




# Recurse into Strawberry source directory
add_subdirectory ( src )
include (CPack)
#link_directories(/usr/local/lib)
#include_directories(include)
#add_library(strawlib common.cpp fasta.cpp contig.cpp read.cpp gff.cpp logger.cpp)
#add_library(estimate.o OBJECT estimate.cpp include/assembly.h)
#target_include_directories(estimate.o PRIVATE)
#add_library(alignments.o OBJECT alignments.cpp include/alignments.h)
#target_include_directories(alignments.o PRIVATE)
#add_library(assembly.o OBJECT assembly.cpp include/assembly.h include/common.h)
#target_include_directories(assembly.o PRIVATE)

#add_executable(Strawberry Strawberry.cpp  $<TARGET_OBJECTS:alignments.o> $<TARGET_OBJECTS:assembly.o> $<TARGET_OBJECTS:estimate.o>) 
#target_link_libraries(Strawberry strawlib bam z pthread )

