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
#include<iostream>
#include "common.h"

class PairedHit;
class ReadHit;
struct CigarOp;

enum Match_t
{
   S_MATCH = 0,
   S_INTRON = 1,
   S_GAP = 2
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
   GenomicFeature(const Match_t& op, uint offset, int len);
   void left(uint left);
   uint left() const;
   void right(uint right);
   uint right() const;

   bool overlap_in_genome(const GenomicFeature& lhs, const GenomicFeature& rhs);

   bool contains(const GenomicFeature& other) const;

   bool properly_contains(const GenomicFeature& other) const;

   static int match_length(const GenomicFeature &op, int left, int right);

   bool operator==(const GenomicFeature & rhs) const;

   bool operator<(const GenomicFeature & rhs) const;

   bool operator!=(const GenomicFeature &rhs) const{
      return !(*this == rhs);
   }

   void printOut(){
      printf("<%d,%d>\n", left(), right());
   }
};

bool readhit_2_genomicFeats(const ReadHit& rh, const std::vector<GenomicFeature> & feats);

/* Reduced representation of GffmRNA using a vector of GenomicFeatures which include intron
 * exon, and unknown*/

class Contig{
   bool _is_ref;
   //std::vector<const PairedHit* > _mates_in_contig;
   std::string _seq;
   RefID _ref_id;
   Strand_t _strand;
   std::string _annotated_trans_id;
   float _mass = 0.0;
public:
   std::vector<GenomicFeature> _genomic_feats;
   Contig() = default;
   Contig(const PairedHit& ph);
   Contig(RefID ref_id, Strand_t strand, const std::vector<GenomicFeature> &feats, bool is_ref);
   const std::string annotated_trans_id() const;
   void annotated_trans_id(std::string str);
   uint left() const;
   uint right() const;
   static bool overlaps_directional(const Contig &lhs, const Contig &rhs);
   bool operator<(const Contig &rhs) const;
   RefID ref_id() const;
   Strand_t strand() const;
   size_t featSize() const;
   float mass() const;
   void mass(float m);
};

#endif /* CONTIG_H_ */
