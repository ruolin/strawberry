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
//#include "qp.h"
#include <chrono>
using namespace std;
void mean_and_sd_insert_size(const vector<int> & vec, double & mean, double &sd){
   double sum = accumulate(vec.begin(), vec.end(), 0.0);
   mean = sum / vec.size();
   cout<<mean<<endl;
   double sq_sum = inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
   sd = std::sqrt(sq_sum / vec.size() - mean * mean);

}

int main(){
   const char *path = "/home/ruolin/Dropbox/Strawberry/Arabidopsis";
   const char *ara_gtf = "/home/ruolin/Dropbox/Strawberry/TAIR10_GFF3_genes-1.gff";
   //const char *human_gtf = "/home/ruolin/Downloads/gencode.v21.annotation.gff3";
   const char *bam_file = "/home/ruolin/Dropbox/Strawberry/accepted_hits.bam";
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
   //QpSolver qps;

   //read_clusters.inspectCluster();
   double mean, sd;
   //const vector<int> & fd = read_clusters._hit_factory->_reads_table._frag_dist;
   //mean_and_sd_insert_size(fd, mean, sd);
   unique_ptr<InsertSize> insert_size(new InsertSize(350, 25));
   read_clusters._insert_size_dist = move(insert_size);
   cout<<"pdf "<<read_clusters._insert_size_dist->truncated_normal_pdf(450)<<endl;
   read_clusters.parseClusters(pFile);


   fclose(pFile);

   auto end = chrono::steady_clock::now();
   auto diff = end - start;
   cout << "Finished in " << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}


