/*
 * contig.h
 *
 *  Created on: Jan 19, 2015
 *      Author: Ruolin Liu
 */

#ifndef CONTIG_H_
#define CONTIG_H_

//#include "GBase.h"
#include<string>
#include<cassert>
#include<vector>
#include<set>
#include<map>
#include<iostream>
#include "common.h"

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
            os<<"exon:["<<gf.left()<<"-"<<gf.right()<<"]\t";
            break;
        case Match_t::S_INTRON:
            os<<"intron:["<<gf.left()<<"-"<<gf.right()<<"]\t";
            break;
        case Match_t::S_GAP:
            os<<"gap:["<<gf.left()<<"-"<<gf.right()<<"]\t";
            break;
        default:
            break;
    }
    return os;
}



bool readhit_2_genomicFeats(const ReadHit& rh, std::vector<GenomicFeature> & feats);

inline auto merge_genomicFeats(const std::vector<GenomicFeature>& feats) {
   std::vector<GenomicFeature> result;

   for(size_t i=0; i<feats.size(); ++i){
      result.push_back(feats[i]);
      GenomicFeature & f = result.back();
      while(i < feats.size() - 1&&
            f._match_op._code == feats[i + 1]._match_op._code)
      {
         if (f._match_op._code == Match_t::S_INTRON) {
            if(f != feats[i + 1]) {
               return std::vector<GenomicFeature>();
            }
         } else {
            if (f.right() < feats[i+1].left()) {
               return std::vector<GenomicFeature>();
            }
            unsigned right = std::max(f.right(), feats[i + 1].right());
            f._match_op._len = right - f.left() + 1;
         }
         ++i;
      }
   }
//   for (auto f: result) {
//      std::cerr<<"output "<<f<<std::endl;
//   }
   return result;
}

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
   bool _is_ref;
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
   ReadID contig_id() const {
      if (_is_ref) return 0;
      else return _contig_id;
   }

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
    os<<"contig<"<<contig.contig_id()<<"> ";
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
