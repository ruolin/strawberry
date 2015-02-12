/*
 * alignment.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: ruolin
 */

#include "alignments.h"
#include "fasta.h"
#ifdef DEBUG
   #include <iostream>
#endif
using namespace std;


bool ClusterFactory::loadRefmRNAs(vector<unique_ptr<GffSeqData>> &gseqs, RefSeqTable &rt,
      const char *seqFile)
{
   //sort gseqs accroding to the observation order in ref_table
   // or if ref_table is empty, initialize it according to gseqs
   if(rt.size() == 0){
      for(uint i=0; i<gseqs.size(); ++i){
         rt.get_id(gseqs[i]->_g_seq_name);
      }
   } else{
      for(uint i = 0; i<gseqs.size(); ++i){
         int idx = gseqs[i]->get_gseq_id();
         if( idx != rt.get_id(gseqs[i]->_g_seq_name) ){
            SMessage("Warning: Sam file and Gff file are not sorted in the same order!\n");
            swap(gseqs[i], gseqs[idx]);
         }
      }
   }
   FaInterface fa_api;
   FaSeqGetter *fsg = NULL;
   if(seqFile != NULL){
      fa_api.initiate(seqFile);
   }
   for(uint i = 0; i<gseqs.size(); ++i){//begin for loop
      GffSeqData * gseq = &(*gseqs[i]);
      int f = 0;
      int r = 0;
      int u = 0;
      RefID ref_id = rt.get_id(gseqs[i]->_g_seq_name);
      GffmRNA *mrna = NULL;
      if(fa_api.hasLoad()){
         delete fsg;
         fsg = NULL;
         fsg = new FaSeqGetter();
         fa_api.load2FaSeqGetter(*fsg,gseqs[i]->_g_seq_name);
         if(fsg == NULL){
            SMessage("Reference sequence %s can not be load!\n",gseqs[i]->_g_seq_name.c_str());
         }
      }
      int f_total = gseqs[i]->_forward_rnas.size();
      int r_total = gseqs[i]->_reverse_rnas.size();
      int u_total = gseqs[i]->_unstranded_rnas.size();
      while(!(f==f_total && r==r_total && u == u_total)){
         char strand = 0;
         if(f <f_total){
            mrna = &(*gseqs[i]->_forward_rnas[f++]);
            strand = '+';
         } else if(r < r_total){
            mrna = &(*gseqs[i]->_reverse_rnas[r++]);
            strand = '-';
         } else{
            mrna = &(*gseqs[i]->_unstranded_rnas[u++]);
            strand = '.';
         }

         vector<GenomicFeature> feats;
         for(uint e = 0; e < mrna->_exons.size(); ++e){
            GffExon& ex = *(mrna->_exons[e]);
            feats.push_back(GenomicFeature(S_MATCH, ex._iv.left(), ex._iv.right()-ex._iv.left()+1));
            if( e + 1 < mrna->_exons.size()){
               GffExon& next_ex  = *(mrna->_exons[e+1]);
               feats.push_back(GenomicFeature(S_INTRON, ex._iv.right()+1, next_ex._iv.left()-1-ex._iv.right() ));
            }
         }
         Contig ref_contig(ref_id, strand, feats, true);
         _ref_mRNAs.push_back(ref_contig);
      }// end while loop
   }//end for loop
   delete fsg;
   fsg = NULL;
   return true;
}

double ClusterFactory::next_valid_alignment(){
   const char* hit_buf=NULL;
   size_t hit_buf_size = 0;
   double raw_mass = 0;
   while(true){
      if(!_hit_factory->nextRecord(hit_buf, hit_buf_size)) break;

      ReadHit tmp;
      if(!_hit_factory->getHitFromBuf(hit_buf, tmp)) continue;
      if(tmp.ref_id() == -1) continue; // unmapped read
      raw_mass += tmp.mass(); // suck in read mass for future if mask_gtf is used.

      if(!_last_hit.left() == 1 || !_last_hit.right() == 0){
         if(_last_hit.ref_id() > tmp.ref_id() ||
               ( _last_hit.ref_id() == tmp.ref_id() && _last_hit.left() > tmp.left())
           )
         {
            const string cur_chr_name = _hit_factory->_ref_table.ref_name(tmp.ref_id());
            const string last_chr_name = _hit_factory->_ref_table.ref_name(_last_hit.ref_id());
            SError("Error:BAM file not sort correctly! The current position is %s:%d and previous position is %s:%d.\n",
                  cur_chr_name.c_str(), tmp.left(), last_chr_name.c_str(), _last_hit.left());
         }
      }

      _last_hit = move(tmp);
      break;
   }
   return raw_mass;
}

double ClusterFactory::rewind_hit(const ReadHit& rh)
{
   double mass = rh.mass();
   _hit_factory->undo_hit();
   return mass;
}

bool ClusterFactory::nextCluster(HitCluster &clusterout)
{
   const ReadHit *bh = NULL;
   while(bh == NULL){
      if(!_hit_factory->recordsRemain()) return false;
      double mass = next_valid_alignment();
   }

}
