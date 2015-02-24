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

enum Op_t {S_MATCH, S_INTRON, S_UNKNOWN};

class GenomicFeature{

public:
   Op_t _code;
   uint _genomic_offset;
   int _genomic_length;
   GenomicFeature(const Op_t& cc, uint offset, int len);
   void left(uint left);
   uint left() const;
   void right(uint right);
   uint right() const;

   static bool overlap_in_genome(const GenomicFeature& lhs, const GenomicFeature& rhs);

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
   std::string _annotated_gene_id;
   std::string _annotated_trans_id;
   RefID _ref_id;
   char _strand;
public:
   std::vector<GenomicFeature> _genomic_feats;
   Contig(RefID ref_id, char strand, const vector<GenomicFeature> &feats, bool is_ref);
   uint left() const;
   uint right() const;
   static bool overlaps_directional(const Contig &lhs, const Contig &rhs);
   bool operator<(const Contig &rhs) const;
   RefID ref_id() const;
   const char strand() const;
};

#endif /* CONTIG_H_ */
