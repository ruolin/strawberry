Strawberry source code current version: 0.8.3 beta. 
==================

What is Strawberry?
==================
Strawberry is a C++ program for fast and accurate ab initio transcript reconstruction and quantification from RNA-seq data. It is written in C++11 and is available as open source software. Strawberry leverages the speed and accuracy of transcript assembly and quantification in such a way that processing 10 million simulated reads (after alignment) requires only 90 seconds using a single thread while achieving over 92% correlation with the ground truth, making it the state-of-the-art method.

Strawberry is a genome-guided transcript-level assembler and quantification tool. It runs on aligned RNA-seq data in BAM format. This allows Strawberry to take advantages of the latest reference genome (if possible, a finished and high-quality one) and stat-of-art splice-awareness aligners. The application of a fast flow network algorithm, for assembly speeds up the construction of transcript-based models. The resulting reduced data representation improves quantification of the different isoforms. Strawberry is also able to account for various sequencing bias that is intrinsic to the RNA-seq experiment. For paired-end RNA-seq, Strawberry empirically infers the insert length distribution from the reads that are mapped uniquely and concordantly. If half-mapped reads exist, Strawberry generates the other ends based on the mapped orientations and the insert length distribution. 

Strawberry consists of two module: assembly module and quantification module. The two modules work in a sequential manner. The assembly module parses aligned reads into splicing graphs, and it uses network flow algorithms to selected the most likely transcripts. The quantification module uses a statistical model, more accurately a Poisson generalized linear model with identity link, to assign ambiguous reads to cpm transcripts. Strawberry simultaneous estimates the transcript abundances and corrects for sequencing bias through the EM steps. 

Documentation
===================
The preprint of Strawberry paper is availble [here](http://biorxiv.org/content/early/2016/03/16/043802).
<hr />

Prerequisites
===================
* A C++11 conformant compiler (currently tested with GCC>=4.7)
* Strawberry uses [CMake](https://cmake.org/)(2.8.12+) build system to check and install dependencies, and to compile and install Strawberry. Make sure the CMake version is at least 2.8.12. 
* [Samtools](http://samtools.sourceforge.net/). Strawberry uses an older version of Samtools(v0.1.19) as a dependent library.   
  * If a different version of Samtools (currently Samtools split into 3 projects now and has been updated to v1.3) has been installed, or your system does not have Samtools installed at all, please let Strawberry automatically download and install Samtools for you. You *DONOT* need to do anything.
  * If you have installed Samtools v0.1.19 and want to use it, please make sure that *libbam.a* is in the Samtools source directory and set an environmental variable SAMTOOLS_ROOT to that location. 
    
    * If you use csh or tcsh, at the shell prompt, enter 
      
     `setenv SAMTOOLS_ROOT /path/to/samtools-source-dir`.
    * If you use bash, at the shell promet, enter 
    
     `export SAMTOOLS_ROOT /path/to/samtools-source-dir`.

Download
========

`git clone https://github.com/ruolin/Strawberry.git`

Precompiled binary is available on the [release page] (https://github.com/ruolin/Strawberry/releases) for the latest version. 

Installation
============
1. `sh cmake.sh`
2. `cd Build`
3. `make`

The executable file will be found in the `Strawberry/bin` directory. 
You can add this directory to your PATH variable to complete the installation. 

Running Strawberry
==================

Running Strawberry is relatively easy. You need to have an alignment file in BAM format. This step can be done using any splice-awareness aligner, e.g. Tophat, GSNAP, STAR. The BAM file needs to be sorted according to the genomic positions. If you use Tophat, the default output is already sorted. For other software, you might have to sort their outputs before running Strawberry. This can be done using Samtools command `samtools sort`.

Type  `strawberry accepted_hits.bam` to run the program on the default parameters. 

For the choice of parameters and their meanings type `strawberry` without any argument for help information. 







