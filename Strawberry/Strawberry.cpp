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
   const char *ara_gtf = "/home/ruolin/Dropbox/Strawberry/TAIR10_GFF3_genes.gff";
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
   read_clusters.ParseClusters();
//   while(true){
//      HitCluster cur;
//      if(read_clusters.nextCluster_refGuide(cur) != -1){
//         if(cur.hasRefmRNAs()){
//            cout<<"number of Ref mRNAs "<<cur._ref_mRNAs.size()<<"\tRef cluster: "\
//                  <<cur.ref_id()<<"\t"<<cur.left()<<"-"<<cur.right()<<"\t"<<cur.size()<<endl;
//            cur.collapseHits();
//            cout<<"number of unique hits\t"<<cur._uniq_hits.size()<<endl;
//         }
//         else{
//            cout<<"Novo cluster number of closed mates: "<<cur.ref_id()<<":"<<
//                  cur.left()<<"-"<<cur.right()<<"\t"<<cur.raw_mass()<<endl;
//         }
//      }
//      else{
//         break;
//      }
//   }

//     vector<int> vec = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2
//     ,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
//     3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4};
//   vector<int> dup;
//   auto last = unique2(vec.begin(), vec.end(), back_inserter(dup));
//   vec.erase(last, vec.end());
//   for(int i: vec)
//      cout<<i<<endl;
   auto end = chrono::steady_clock::now();
   auto diff = end - start;
   cout << "Finished in " << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}

//72339-74096
//73931-74737


