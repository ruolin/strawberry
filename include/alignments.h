/*
 * alignment.h
 *
 *  Created on: Nov 5, 2014
 *      Author: ruolin
 */

#ifndef STRAWB_ALIGNMENTS_H_
#define STRAWB_ALIGNMENTS_H_
#include <list>
#include <map>
#include <climits>
#include <atomic>
#include "read.h"
#include "contig.h"
#include "gff.h"
class Sample;
class FaInterface;
class FaSeqGetter;
using IntronMap = std::map<std::pair<uint,uint>,IntronTable>;

class HitCluster{
   friend Sample;
   uint _leftmost;
   uint _rightmost;
   int _plus_strand_num_hits;
   int _minus_strand_num_hits;
   Strand_t _first_encounter_strand;
   int _id;
   RefID _ref_id;
   bool _final; // HitCluster is finished
   double _raw_mass;
   Strand_t _strand;
   std::unordered_map<ReadID, list<PairedHit>> _open_mates;
   std::vector<PairedHit> _hits;
   std::vector<PairedHit> _uniq_hits;
   std::vector<Contig*> _ref_mRNAs; // the actually objects are owned by Sample
   std::vector<GenomicFeature> _introns;
   std::vector<float> _dep_of_cov;
   //std::map<std::pair<int,int>,int> _current_intron_counter;
public:
   //static const int _kMaxGeneLen = 1000000;
   //static const int _kMaxFragPerCluster = 100000;
   static const int _kMinFold4BothStrand = 10;
   HitCluster();
   RefID ref_id() const;
   void ref_id(RefID id);
   uint left() const;
   void left(uint left);
   void right(uint right);
   uint right() const;
   int size() const;
   int len() const;
   Strand_t ref_strand() const;
   Strand_t guessStrand() const;
   Strand_t strand() const{
      return _strand;
   }
   bool addHit(const PairedHit &hit);
   void count_current_intron(const std::vector<std::pair<uint,uint>>& cur_intron);
   void setBoundaries();
   void clearOpenMates();
   bool addOpenHit(ReadHitPtr hit, bool extend_by_hit, bool extend_by_partner);
   int collapseHits();
   bool overlaps(const HitCluster& rhs) const;
   bool hasRefmRNAs() const {
      return _ref_mRNAs.size() > 0;
   }
   void addRefContig(Contig *contig);
   int numOpenMates() const{
      return _open_mates.size();
   }
   void addRawMass(double m){
      _raw_mass += m;
   }
   void subRawMass(double m){
      _raw_mass -= m;
   }
   double raw_mass() const{
      return _raw_mass;
   }
   double collapse_mass() const;
   //void count_current_intron(const ReadHit & hit);
   //bool current_intron_is_reliable() const;
   bool see_both_strands();

   void reset(uint old_left,
              uint old_right,
              RefID old_ref_id,
              int old_plus_hit,
              int old_minus_hit,
              Strand_t old_strand,
              double old_mass
              );

   void reset(uint old_left,
              uint old_right,
              RefID old_ref_id
              );

};

class Sample{
   bool _is_inspecting;
   int _num_cluster = 0;
   //uint _prev_pos = 0;
   RefID _prev_hit_ref_id = -1; //used to judge if sam/bam is sorted.
   uint _prev_hit_pos = 0; //used to judge if sam/bam is sorted.
   size_t _refmRNA_offset;
   bool _has_load_all_refs;
   string _current_chrom;
   int _total_mapped_reads = 0;
   double compute_doc(const uint left,
                     const uint right,
                     const vector<Contig> & hits,
                     vector<float> &exon_doc,
                     IntronMap &intron_doc,
                     uint smallOverhang);
public:
   shared_ptr<HitFactory> _hit_factory;
   shared_ptr<InsertSize> _insert_size_dist =nullptr;
   shared_ptr<FaSeqGetter> _fasta_getter = nullptr;
   shared_ptr<FaInterface> _fasta_interface = nullptr;

   std::vector<Contig> _ref_mRNAs; // sort by seq_id in reference_table
   Sample(shared_ptr<HitFactory> hit_fac):
      _refmRNA_offset(0),
      _has_load_all_refs(false),
      _hit_factory(move(hit_fac))
   {}
   //int max_inner_dist() const;
   int total_mapped_reads() const;
   bool load_chrom_fasta(RefID seq_id);
   bool hasLoadRefmRNAs() const {
      return _ref_mRNAs.size() > 0;
   }
   bool loadRefFasta(RefSeqTable &rt, const char *seqFile = NULL);
   bool loadRefmRNAs(vector<unique_ptr<GffSeqData>> &gseqs, RefSeqTable &rt);
   int addRef2Cluster(HitCluster &clusterOut);
   void reset_refmRNAs();
   double next_valid_alignment(ReadHit& readin);
   double rewindHit(const ReadHit& rh);
   int nextCluster_denovo(HitCluster &clusterOut,
                           uint next_ref_start_pos = UINT_MAX,
                           RefID next_ref_start_ref=INT_MAX);

   int nextCluster_refGuide(HitCluster & clusterOut);
   void rewindReference(HitCluster &clusterOut, int num_regress);

   void mergeClusters(HitCluster & dest, HitCluster &resource);
   //void compute_doc_4_cluster(const HitCluster & hit_cluster, vector<float> &exon_doc,
                              //map<pair<uint,uint>,IntronTable>& intron_counter, uint &small_overhang);

   void filter_intron(uint cluster_left, vector<float> &exon_doc, IntronMap& intron_counter);
   void procSample(FILE *f, FILE *log);
   void inspectSample(FILE *log);
   void finalizeAndAssemble(const RefSeqTable & ref_t, shared_ptr<HitCluster> cluster, FILE *f, FILE *plogfile);
};



bool hit_lt_cluster(const ReadHit& hit, const HitCluster& cluster, uint olap_radius);
bool hit_gt_cluster(const ReadHit& hit, const HitCluster& cluster, uint olap_radius);
bool hit_complete_within_cluster(const PairedHit& hit, const HitCluster& cluster, uint olap_radius);
#endif /* STRAWB_ALIGNMENTS_H_ */
