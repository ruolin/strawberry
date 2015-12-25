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
#include<iostream>
#include "common.h"

class PairedHit;
class ReadHit;
class RefSeqTable;
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
   uint _genomic_offset;
   MatchOp _match_op;
   double _avg_cov; // average coverage
   GenomicFeature(const Match_t& op, uint offset, int len);
   int len() const;
   void left(uint left);
   uint left() const;
   void right(uint right);
   uint right() const;
   void avg_doc(double coverage);
   double avg_doc() const;

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

   void printOut(){
      printf("<%d,%d>\n", left(), right());
   }
   static void mergeFeatures(const vector<GenomicFeature> & feats, vector<GenomicFeature> & result);

};


bool readhit_2_genomicFeats(const ReadHit& rh, const std::vector<GenomicFeature> & feats);

/* Reduced representation of GffmRNA using a vector of GenomicFeatures which include intron
 * exon, and unknown*/

class Contig{
   //std::vector<const PairedHit* > _mates_in_contig;
   RefID _ref_id;
   std::string _seq;
   Strand_t _strand;
   std::string _annotated_trans_id;
   double _mass = 0.0;
   SingleOrit_t _single_read_orit = SingleOrit_t::NotSingle;
public:
   std::vector<GenomicFeature> _genomic_feats;
   bool _is_ref;
   Contig() = default;
   Contig(const PairedHit& ph);
   Contig(RefID ref_id, Strand_t strand, const std::vector<GenomicFeature> &feats, double mass,bool is_ref);
   //Contig(const ExonBin& eb);
   const std::string annotated_trans_id() const;
   void annotated_trans_id(std::string str);
   uint left() const;
   uint right() const;
   uint gap_left() const; // left coordinate of gap if exists; otherwise return 0
   uint gap_right() const; // left coordiante of gap if exists; otherwise return 0
   int exonic_length() const;
   static int exonic_overlaps_len(const Contig &iso,
         const uint left,
         const uint right);
   static bool overlaps_directional(const Contig &lhs, const Contig &rhs);
   static bool is_contained_in(const Contig &small, const Contig &large);
   static bool is_compatible(const Contig &read, const Contig &isoform);
   //static int infer_inner_dist(const Contig &isoform, const Contig &hit);

   bool operator<(const Contig &rhs) const;
   RefID ref_id() const;
   Strand_t strand() const;
   size_t featSize() const;
   bool is_single_read() const;
   float mass() const;
   void mass(float m);
   void print2gtf(FILE *pFile,
                  const RefSeqTable &ref_lookup,
                  const string fpkm,
                  const string tpm,
                  int gene_id, int tscp_id);

   SingleOrit_t single_read_orit() const;
};

#endif /* CONTIG_H_ */
