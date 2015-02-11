/*
 * alignment.h
 *
 *  Created on: Nov 5, 2014
 *      Author: ruolin
 */

#ifndef STRAWB_ALIGNMENTS_H_
#define STRAWB_ALIGNMENTS_H_
#include "read.h"
#include "contig.h"
#include "gff.h"
class HitCluster{
   int _leftmost = INT_MAX;
   int _rightmost = -1;
   int _id;
   RefID _ref_id = 0;
   bool _final = false;
   double _raw_mass = 0.0;
   static int _next_id;
public:
   std::vector<PairedHit> _hits;
   std::vector<PairedHit> _non_redundant;
   std::vector<Contig*> _ref_contigs; // the actually objects are owned by ClusterFactory
   std::vector<GenomicFeature> _introns;
   HitCluster(): _id(++_next_id){}
};


#endif /* STRAWB_ALIGNMENTS_H_ */


class ClusterFactory{
   unique_ptr<HitFactory> _hit_factory;
   int _num_cluster = 0;
   std::vector<Contig>::iterator _next_ref_mRNA;
   uint _prev_pos = 0;
   RefID _prev_ref_id = 0;
public:
   ReadHit _last_hit; // Record the hit that next_valid_alignment() read in.
   std::vector<Contig> _ref_mRNAs;
   ClusterFactory(unique_ptr<HitFactory> hit_fac):
      _hit_factory(move(hit_fac))
   {}
   bool nextCluster(HitCluster & clusterOut);
   bool loadRefmRNAs(vector<unique_ptr<GffSeqData>> &gseqs, RefSeqTable &rt, const char *seqFile = NULL);
   int num_refmRNA() const {
      return _ref_mRNAs.size();
   }
   double next_valid_alignment();
   double rewind_hit(const ReadHit& rh);
};
