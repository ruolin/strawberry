cmake_minimum_required(VERSION 2.8.12)

project( strawberry LANGUAGES CXX)
set( CMAKE_CXX_STANDARD 14 )
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

set(CPACK_PACKAGE_VERSION "1.0.5")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "5")
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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -funroll-loops -fPIC -fomit-frame-pointer -Wall -Wno-reorder -Wno-unused-variable -Wreturn-type -Werror=return-type")
set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Ofast")
#set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -g")
add_definitions(-DELPP_THREAD_SAFE)

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
                            "g++ >= 5.3.0")
    elseif (NOT (GCC_VERSION VERSION_GREATER 4.8.2 OR GCC_VERSION VERSION_EQUAL 4.8.2))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 5.3.0 or greater.")
    endif ()

    set (GCC TRUE)

    set (PTHREAD_LIB "pthread")

    # If we're on Linux (i.e. not OSX) and we're using
    # gcc, then set the -static-libstdc++ flag
    if (NOT APPLE)
        set( CMAKE_EXE_LINKER_FLAGS " -static-libgcc -static-libstdc++ -L. ")
        #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")
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
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)

if(COMPILER_SUPPORTS_CXX14)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
#elseif(COMPILER_SUPPORTS_CXX11)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#elseif(COMPILER_SUPPORTS_CXX0X)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
   message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a C++ compiler that supports C++14.")
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
ExternalProject_Add(eigen
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
  URL http://bitbucket.org/eigen/eigen/get/3.3.3.tar.bz2 
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen3
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
)


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
  BUILD_COMMAND make
  INSTALL_COMMAND ""
)

set(BAM_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools-0.1.19/libbam.a)
endif()




# Recurse into Strawberry source directory
add_subdirectory ( src )
include (CPack)

