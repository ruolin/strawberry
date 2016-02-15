What is Strawberry?
==================
Strawberry is a C++ program for fast yet accurate ab initio transcript reconstruction and quantification from RNA-seq data. It is written in C++11, and is available as open source software. It leverages the speed and the accuracy in such a way that it finishes assembling and quantifing 10 million aligned artifical reads within 2 minutes on a 16GB RAM and core-i7 desktop while achieving over 90% correlation with the ground truth. This accuracy is very close to a very differnt yet state-of-art known transcript quantification approach, RSEM. 

Strawberry is a genome-guided transcript-level assembler and quantification tool. It runs on aligned RNA-seq data in BAM format. Strawberry is specially designed for paired-end data based on our belief that it is advantageous to use paired-end RNA-seq over single-end data. It empirically infers the insert length distribution from the reads that are mapped uniquely and concordantly. For those halfmapped reads, it simulates their other end based on the mapped orientation and the insert length distribution. 

Strawberry consists of two module: assembly module and quantifcation module. The two modules work in a sequential manner. The assembly module parses aligned reads into splicing graphs and it uses network flow algorithms to selected the most possible transcripts. The quantification module uses a statistcal model, more accuratly a Poisson generalized linear model with identity link, to assign ambiguous reads to possible transcripts. Strawberry simultaneous estimates the transcript abundances and corrects for sequencing bias through the EM steps. 
<hr />

Prerequisites
===================
* A C++11 conformant compiler (currently tested with GCC>=4.7)
* Strawberry uses the [CMake](https://cmake.org/)(2.8.12+) to build system to check, fetch and install dependencies, and to compile and install Strawberry. Make sure the Cmake version is at least 2.8.12. 
* [Samtools](http://samtools.sourceforge.net/). Strawberry uses an older version of Samtools(v0.1.19) as a dependent library.   
  * If your system installed a different version (currently Samtools split into 3 projects now and has been updated to v1.3) or your system does not have Samtools installed at all, please let Strawberry automatically download and install Samtools for you. You DONOT need to do anything
  * If you have install Samtools v0.1.19 and want to use it, please make sure that *libbam.a* is in the Samtools source directory and set an environmental variable SAMTOOLS_ROOT to that location. 
    
    * If you use csh or tcsh, at the shell prompt, enter 
      
     `setenv SAMTOOLS_ROOT /path/to/samtools-source-dir`
    * If you use bash, at the shell promet, enter 
    
     `export SAMTOOLS_ROOT /path/to/samtools-source-dir`

Download
========

`git clone https://github.com/ruolin/Strawberry.git`

Installation
============
1. `sh cmake.sh`
2. `cd Build`
3. `make`

The executable file will be found in the `Strawberry/bin` directory. 
You can add this directory to your PATH variable to complete the installation. 

Running Strawberry
==================

Running Strawberry is fairly easy. You need to have an alignemnt file in BAM format. This step can be done using any splice-awareness aligner, e.g. Tophat, GSNAP, STAR. The BAM file has to be sorted according to the genomic positions. If you use Tophat, the default output is already sorted. For other software, you might have to sort their outputs before running Strawberry. This can be done using Samtools command `samtools sort`

type  `strawberry accepted_hits.bam` to run the program on the default parameters. 

For the choice of parameters and their meanings type `strawberry` without any argument for help informations. 







