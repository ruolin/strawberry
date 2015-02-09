/*
 * alignment.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: ruolin
 */

#include "alignments.h"
#include "fasta.h"
using namespace std;


bool ClusterFactory::loadRefmRNA(const vector<unique_ptr<GffSeqData>> &gseqs, RefSeqTable &rt,
      const char *seqFile)
{
   FaInterface fa_api(seqFile);
   FaSeqGetter fsg;
   fa_api.load2FaSeqGetter(fsg, "mitochondria");

   //sort gseqs accroding to the observation order in ref_table
   // or if ref_table is empty, initialize it according to gseqs
   if(rt.size() == 0){
      for(int i=0; i<gseqs.size(); ++i)
         rt.get_id(gseqs[i]->_g_seq_name);
   } else{
      for(int i = 0; i<gseqs.size(); ++i){
         int idx = gseqs[i]->get_gseq_id();
         if( idx != rt.get_id(gseqs[i]->_g_seq_name) ){
            SMessage("Warning: Sam file and Gff file are not sorted in the same order!\n");
            swap(gseqs[i], gseqs[idx]);
         }
      }
   }

   for(int i = 0; i<gseqs.size(); ++i){
      GffSeqData * gseq = &(*gseqs[i]);
      int f = 0;
      int r = 0;
      int u = 0;

   }

}
