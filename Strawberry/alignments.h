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
   uint _leftmost = UINT_MAX;
   uint _rightmost = 0;
   int _id;
   RefID _ref_id = -1;
   bool _final = false; // HitCluster is finished
   double _raw_mass = 0.0;
   static int _next_id ;
   std::unordered_map<ReadID, unique_ptr<PairedHit>> _open_mates;
public:
   static const int _kMaxGeneLen = 1000000;
   static const int _kMaxFragPerCluster = 100000;
   std::vector<PairedHit> _hits;
   std::vector<PairedHit> _non_redundant;
   std::vector<Contig*> _ref_mRNAs; // the actually objects are owned by ClusterFactory
   std::vector<GenomicFeature> _introns;
   HitCluster();
   RefID ref_id() const;
   void ref_id(RefID id);
   uint left() const;
   void left(uint left);
   void right(uint right);
   uint right() const;
   int size() const;
   bool addHit(const PairedHit &hit);
   bool addOpenHit(unique_ptr<ReadHit> hit, bool extend_by_hit);
   bool hasRefmRNAs() const {
      return _ref_mRNAs.size() > 0;
   }
   void addRefContig(Contig *contig);
   int numOpenMates() const{
      return _open_mates.size();
   }
};


#endif /* STRAWB_ALIGNMENTS_H_ */


class ClusterFactory{
   static const int _kMinOlapDist = 50;
   unique_ptr<HitFactory> _hit_factory;
   int _num_cluster = 0;
   uint _prev_pos = 0;
   RefID _prev_hit_ref_id = -1; //used to judge if sam/bam is sorted.
   uint _prev_hit_pos = 0; //used to judge if sam/bam is sorted.
   int _refmRNA_offset;
   bool _has_load_all_refs;
public:
   std::vector<Contig> _ref_mRNAs; // sort by seq_id in reference_table
   ClusterFactory(unique_ptr<HitFactory> hit_fac):
      _hit_factory(move(hit_fac)),
      _refmRNA_offset(0),
      _has_load_all_refs(false)
   {}

   bool loadRefmRNAs(vector<unique_ptr<GffSeqData>> &gseqs, RefSeqTable &rt, const char *seqFile = NULL);
   bool hasLoadRefmRNAs() const {
      return _ref_mRNAs.size() > 0;
   }
   double next_valid_alignment(ReadHit& readin);
   double rewindHit(const ReadHit& rh);
   int addRef2Cluster(HitCluster &clusterOut);
   int nextCluster_denovo(HitCluster &clusterOut,
                           uint next_ref_start_pos = UINT_MAX,
                           RefID next_ref_start_ref=INT_MAX);

   int nextCluster_refGuide(HitCluster & clusterOut);
   void rewindReference(HitCluster &clusterOut, int num_regress);
};

bool hit_lt_cluster(const ReadHit& hit, const HitCluster& cluster);
bool hit_gt_cluster(const ReadHit& hit, const HitCluster& cluster);
