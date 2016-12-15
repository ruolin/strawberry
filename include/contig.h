/*
 * contig.h
 *
 *  Created on: Jan 19, 2015
 *      Author: ruolin
 */

#ifndef CONTIG_H_
#define CONTIG_H_

//#include "GBase.h"
#include<string>
#include<cassert>
#include<vector>
#include<set>
#include<map>
#include<algorithm>
#include<iostream>
#include "common.h"
#include "fasta.h"


class PairedHit;
class ReadHit;
class RefSeqTable;
class Contig;
struct CigarOp;

enum Match_t
{
   S_MATCH = 0,
   S_INTRON = 1,
   S_GAP = 2
};

enum SingleOrit_t
{
   Forward = 0,
   Reverse = 1,
   NotSingle = 2
};

struct MatchOp{
   Match_t _code : 2;
   uint32_t _len : 30;
   MatchOp(Match_t code, uint32_t len): _code(code), _len(len){}
   bool operator==(const MatchOp & rhs) const{
      return _code == rhs._code && _len == rhs._len;
   }
};

class GenomicFeature{

public:
    using TDepth = int;
   uint _genomic_offset;
   MatchOp _match_op;
   double _avg_cov; // average coverage
   GenomicFeature(const Match_t& op, uint offset, int len);
   GenomicFeature(): GenomicFeature(Match_t::S_MATCH, 0u, 0) {};
   int len() const;
   void left(uint left);
   uint left() const;
   void right(uint right);
   uint right() const;
   int depth() const {return 1;}
   void avg_doc(double coverage);
   double avg_doc() const;
   bool compatible_2_read(const Contig& read) const;
   static bool overlaps(const GenomicFeature& lhs, const GenomicFeature& rhs);
   static int overlap_len(const GenomicFeature &lhs, const GenomicFeature &rhs);
   static int overlap_len_in_genome(const GenomicFeature& feat, const uint left, const uint right);
   static bool overlap_in_genome(const GenomicFeature & feat, const uint left, const uint right);
   bool contains(const GenomicFeature& other, int small_extent = 0) const;

   bool properly_contains(const GenomicFeature& other) const;

   static int match_length(const GenomicFeature &op, int left, int right);

   friend bool operator==(const GenomicFeature &lhs, const GenomicFeature & rhs);

   bool operator<(const GenomicFeature & rhs) const;

   friend bool operator!=(const GenomicFeature &lhs, const GenomicFeature &rhs);
    friend std::ostream& operator<<(std::ostream&, const GenomicFeature& );

   void printOut(){
      printf("<%d,%d>\n", left(), right());
   }
   static void mergeFeatures(const std::vector<GenomicFeature> & feats, std::vector<GenomicFeature> & result);
};

inline std::ostream& operator<<(std::ostream& os, const GenomicFeature& gf){
    switch(gf._match_op._code){
        case Match_t::S_MATCH:
            os<<"exon:["<<gf.left()<<"-"<<gf.right()<<"]("<<gf.len()<<")\t";
            break;
        case Match_t::S_INTRON:
            os<<"intron:["<<gf.left()<<"-"<<gf.right()<<"]("<<gf.len()<<")\t";
            break;
        case Match_t::S_GAP:
            os<<"gap:["<<gf.left()<<"-"<<gf.right()<<"]("<<gf.len()<<")\t";
            break;
        default:
            break;
    }
    return os;
}



bool readhit_2_genomicFeats(const ReadHit& rh, const std::vector<GenomicFeature> & feats);

/* Reduced representation of GffmRNA using a vector of GenomicFeatures which include intron
 * exon, and unknown*/

class Contig{
   //std::vector<const PairedHit* > _mates_in_contig;
   RefID _ref_id;
   ReadID _contig_id ;

   //std::string _seq;
   Strand_t _strand;
   std::string _annotated_trans_id;
   std::string _parent_id;
   double _mass = 0.0;
   SingleOrit_t _single_read_orit = SingleOrit_t::NotSingle;
public:
   bool _is_ref = false;
   std::vector<GenomicFeature> _genomic_feats;
   Contig() = default;

   Contig(const PairedHit& ph);

   Contig(
           RefID ref_id,
           ReadID cid,
           Strand_t strand,
           double mass,
           std::vector<GenomicFeature> feats,
           bool is_ref):
           _ref_id(ref_id),
           _contig_id(cid),
           _strand(strand),
           _mass(mass),
           _genomic_feats(feats),
           _is_ref(is_ref)
   {
      assert(_genomic_feats.front()._match_op._code == Match_t::S_MATCH);
      assert(_genomic_feats.back()._match_op._code == Match_t::S_MATCH);
   }

   const std::string annotated_trans_id() const;
   void annotated_trans_id(std::string str);

   decltype(auto) parent_id() {return (_parent_id);}
   decltype(auto) parent_id() const {return (_parent_id);}
   uint left() const;
   uint right() const;
   uint gap_left() const; // left coordinate of gap if exists; otherwise return 0
   uint gap_right() const; // left coordiante of gap if exists; otherwise return 0
   int exonic_length() const;
   decltype(auto) id() { return (_contig_id);}
   decltype(auto) id() const { return (_contig_id);}

   static int exonic_overlaps_len(const Contig &iso,
         const uint left,
         const uint right);
   static bool overlaps_directional(const Contig &lhs, const Contig &rhs);
   static bool is_contained_in(const Contig &small, const Contig &large);
   static bool is_compatible(const Contig &read, const Contig &isoform);
   static bool is_compatible(const Contig &isoform, const GenomicFeature& feat);
   //static int infer_inner_dist(const Contig &isoform, const Contig &hit);

   static uint read_start_from_iso(const Contig &iso, const Contig& hit);
   static int fragment_len(const Contig& read, const Contig& iso);

   static int relative_pos_from_left(const Contig& read, const Contig& iso){
      if (is_compatible(read, iso)) {
         int start = 0;
         assert(read.left() >= iso.left());
         for (const auto& gfeat : iso._genomic_feats) {
            if (gfeat._match_op._code == Match_t::S_MATCH) {
               if (read.left() >= gfeat.left()) {
                  if (read.left() > gfeat.right()) start += gfeat.right() - gfeat.left();
                  else start += read.left() - gfeat.left() + 1;
               }
            }
         }
         return start;
      }
      else {
         return 0;
      }
   }

   static int relative_pos_from_right(const Contig& read, const Contig& iso){
      if (is_compatible(read, iso)) {
         int start = 0;
         assert(read.right() <= iso.right());
         for (auto rit = iso._genomic_feats.crbegin(); rit != iso._genomic_feats.crend(); ++rit) {
            if (rit->_match_op._code == Match_t::S_MATCH) {
               if (read.right() <= rit->right()) {
                  if (read.right() < rit->left()) start += rit->right() - rit->left() +1;
                  else start += rit->right() - read.right() + 1;
               }
            }
         }
         return start;
      }
      else {
         return 0;
      }
   }

   static float contig_gc_content(const Contig& contig, const std::shared_ptr<FaSeqGetter> &fa_getter) {
      return contig_gc_content(contig._genomic_feats, fa_getter);
   }


   static float contig_gc_content(const std::vector<GenomicFeature>& gfeats, const std::shared_ptr<FaSeqGetter> &fa_getter) {
      assert(fa_getter != nullptr);
      std::vector<int> lengths;
      std::vector<float> gcs;
      for (const auto& gf : gfeats) {
         if (gf._match_op._code == Match_t::S_MATCH) {
            float gc = gc_content(fa_getter->fetchSeq(gf.left(), gf.len()));
            gcs.push_back(gc);
            lengths.push_back(gf.len());
         }
      }
      float sum = 0.0;
      for (size_t i = 0; i < lengths.size(); ++i) {
         sum += gcs[i] * lengths[i];
      }
      return sum / std::accumulate(lengths.begin(), lengths.end(), 0);
   }

   static std::vector<double> start_site_dist(const Contig & iso, const std::vector<Contig> & hits);
   bool operator<(const Contig &rhs) const;
   friend bool operator==(const Contig &lhs, const Contig & rhs);

   RefID ref_id() const;
   Strand_t strand() const;
   size_t featSize() const;
   bool is_single_read() const;
   float mass() const;
   void mass(float m);
   void print2gtf(FILE *pFile,
                  const RefSeqTable &ref_lookup,
                  const std::string fpkm,
                  const std::string tpm,
                  std::string gene_id, std::string tscp_id) const;

   SingleOrit_t single_read_orit() const;
   double avg_doc() const;
};

inline std::ostream& operator<<(std::ostream& os, const Contig& contig){
    os<<"contig<"<<contig.id()<<"> ";
    os<<contig.ref_id()<<":"<<contig.left()<<"-"<<contig.right()<<"\t";
    for (const auto& gf: contig._genomic_feats) {
        os<<gf;
    }
    os<<std::endl;
    return os;
}
//template<typename TContig>
//class ContigGroup<TContig>{
//private:
//    vector<TContig> contigs_;
//public:
//    decltype(auto) contigs() {return contigs_;}
//};

#endif /* CONTIG_H_ */
