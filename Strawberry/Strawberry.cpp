//============================================================================
// Name        : Strawberry.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : For Strawberry
//============================================================================

#include <iostream>
#include <algorithm>
#include "fasta.h"
#include "gff.h"
#include "alignments.h"
#include "logger.hpp"
#include <chrono>
using namespace std;


int main(){
   const char *path = "/home/ruolin/Dropbox/Strawberry/Arabidopsis";
   const char *ara_gtf = "/home/ruolin/Dropbox/Strawberry/TAIR10_GFF3_genes-1.gff";
   //const char *human_gtf = "/home/ruolin/Downloads/gencode.v21.annotation.gff3";
   const char *bam_file = "/home/ruolin/Dropbox/Strawberry/RD25.high.diffMean_r1.concordant_uniq.sort.bam";
   //const char *bam_file = "/home/ruolin/Dropbox/Strawberry/WetFT1.sm.bam";

   //FaInterface fa_api(path);
   //FaSeqGetter fsg;
   //fa_api.load2FaSeqGetter(fsg, "mitochondria");
   //cout<<"success\t"<<fsg.loadSeq()<<endl;
   //cout<<fsg.fetchSeq(80,4)<<endl;
   auto start = chrono::steady_clock::now();
   GffReader greader(ara_gtf);
   greader.readAll();
   greader.reverseExonOrderInMinusStrand();
   ReadTable read_table;
   RefSeqTable ref_seq_table(true);
   unique_ptr<HitFactory> hf(new BAMHitFactory(bam_file, read_table, ref_seq_table));
   hf->inspect_header();
   ClusterFactory read_clusters(move(hf));
   read_clusters.loadRefmRNAs(greader._g_seqs, ref_seq_table, path);
   FILE *pFile;
   pFile = fopen("assembled_transcripts.gtf", "w");
   read_clusters.ParseClusters(pFile);
   fclose(pFile);

   auto end = chrono::steady_clock::now();
   auto diff = end - start;
   cout << "Finished in " << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}

//72339-74096
//73931-74737


