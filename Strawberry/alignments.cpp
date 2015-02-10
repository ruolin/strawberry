/*
 * alignment.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: ruolin
 */

#include "alignments.h"
#include "fasta.h"
using namespace std;


bool ClusterFactory::loadRefmRNA(vector<unique_ptr<GffSeqData>> &gseqs, RefSeqTable &rt,
      const char *seqFile)
{
   //sort gseqs accroding to the observation order in ref_table
   // or if ref_table is empty, initialize it according to gseqs
   if(rt.size() == 0){
      for(uint i=0; i<gseqs.size(); ++i)
         rt.get_id(gseqs[i]->_g_seq_name);
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
      }// end while loop
   }//end for loop
   delete fsg;
   fsg = NULL;
}
