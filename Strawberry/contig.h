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
#include "read.h"

enum Match_t
{
   S_MATCH,
   S_INTRON,
   S_UNKNOWN
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

   void printOut(){
      printf("<%d,%d>\n", left(), right());
   }
};



/* Reduced representation of GffmRNA using a vector of GenomicFeatures which include intron
 * exon, and unknown*/

class Contig{
   bool _is_ref;
   std::vector<const PairedHit* > _mates_in_contig;
   std::string _seq;
   RefID _ref_id;
   Strand_t _strand;
   std::string _annotated_trans_id;
   std::vector<GenomicFeature> _genomic_feats;
public:
   const std::string annotated_trans_id() const;
   void annotated_trans_id(std::string str);
   Contig(RefID ref_id, Strand_t strand, const vector<GenomicFeature> &feats, bool is_ref);
   uint left() const;
   uint right() const;
   static bool overlaps_directional(const Contig &lhs, const Contig &rhs);
   bool operator<(const Contig &rhs) const;
   RefID ref_id() const;
   Strand_t strand() const;
   size_t featSize() const;
};

#endif /* CONTIG_H_ */
