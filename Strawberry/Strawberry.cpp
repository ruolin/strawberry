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
#include <chrono>
using namespace std;



int main(){
   const char *path = "/home/ruolin/Dropbox/Arabidopsis";
   const char *ara_gtf = "/home/ruolin/Dropbox/Strawberry/TAIR10_GFF3_genes-1.gff";
   //const char *human_gtf = "/home/ruolin/Downloads/gencode.v21.annotation.gff3";
   const char *bam_file = "/home/ruolin/Dropbox/StrawberryWetFT1.sm.bam";
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
   cout<<ref_seq_table.size()<<endl;
   ClusterFactory read_clusters(move(hf));
   read_clusters.loadRefmRNAs(greader._g_seqs, ref_seq_table, path);
   while(true){
      double mass =read_clusters.next_valid_alignment();
      if(mass == 0.0){
         break;
      }
   };
   auto end = chrono::steady_clock::now();
   auto diff = end - start;
   cout << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}


/*

int main() {

   //GenomicInterval iv3("chr1", 0, 120);
   //GffGene g1(iv3, "gene1");
   //foo(g1);
   //cout<<g1._transcripts["transcript2"]->_trans_id<<endl;
   FILE* fh;
   GffReader *gffr = new GffReader("/home/ruolin/Dropbox/Strawberry/test.gff");
   gffr->readAll();


   BAMHitFactory * bh = new BAMHitFactory("/home/ruolin/Dropbox/Strawberry/CoolL1T1.sm.bam");


   while(true){

      const char* hit_buf;
      ReadHit tmp;
      size_t hit_buf_size = 0;
      bh->nextRecord(hit_buf, hit_buf_size);
      bh->getHitFromBuf(hit_buf, tmp);
      GenomicInterval hit_iv(tmp._ref_id.c_str(), tmp._start, tmp.getEnd(), tmp._strand);
//      cout<<hit_iv.getChrom()<<"\t"<<hit_iv.getStart()<<"\t"<<hit_iv.getEnd()<<"\t"<<hit_iv.getStrand()<<endl;
      for (auto gene : gffr->_genes){
         if(gene.second->_transcript_num == 1 && gene.second->_gene_iv.getStrand() == strand_plus){
            for(auto transcript : gene.second->_transcripts){
               transPtr now = transcript.second;
               if (transcript.second->_trans_iv.contain(&hit_iv)){
                  for(auto i=0; i< transcript.second->_exons.size();++i){
                     if (transcript.second->_exons[i]->_exon_iv.overlap(&hit_iv)){
                        if (hit_iv.getStrand() == '+' ){
                           int read_start = hit_iv.getEnd() - transcript.second->_exons[i]->_exon_iv.getStart();
                           for ( auto j=0; j<i; ++j)
                              read_start += transcript.second->_exons[j]->_exon_iv.len();
                           cout<< now->_trans_id<<"\t"<< tmp._read_id<<"\t"<<i+1<<"\t"<<read_start<<endl;
                           break;
                        }
                     } // end exon loop
                  }
               }
            } // end transcript loop
         }
      } // end gene loop
      if(tmp._start > 300000) break;
   }


   return EXIT_SUCCESS;

}

*/


/*  Sample code for reading bam file;
 *    samfile_t* _hit_file;
   _hit_file = samopen("WetFT1.sm.bam", "rb", 0);
   bam1_t _next_hit;
   int start;
   int end;
   int target_id;
   while(true)
   {
      int bytes_read = samread(_hit_file, &_next_hit);
      int a = bam1_cigar(&_next_hit)[0] >> BAM_CIGAR_SHIFT;
      if(a<75){
         start = _next_hit.core.pos;
         end = _next_hit.core.pos + _next_hit.core.l_qseq;
         int b = bam1_cigar(&_next_hit)[0] & BAM_CIGAR_MASK;
         target_id = _next_hit.core.tid;
         break;
      }
   //const char*& buf;
   //buf = (const char*)&_next_hit
   }
   string text_name = _hit_file->header->target_name[target_id];
 *
 */
