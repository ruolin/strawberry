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
   void g_left(uint left);
   uint g_left() const;
   void g_right(uint right);
   uint g_right() const;

   static bool overlap_in_genome(const GenomicFeature& lhs, const GenomicFeature& rhs);

   bool contains(const GenomicFeature& other) const;

   bool properly_contains(const GenomicFeature& other) const;

   static int match_length(const GenomicFeature &op, int left, int right);

   bool operator==(const GenomicFeature & rhs) const;

   bool operator<(const GenomicFeature & rhs) const;

   void printOut(){
      printf("<%d,%d>\n", g_left(),g_right());
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
   Contig(RefID ref_id, char strand, vector<GenomicFeature> &feats, bool is_ref);
   int left() const;
   int right() const;
   bool operator<(const Contig &rhs) const;
   RefID get_ref_id() const;
};

#endif /* CONTIG_H_ */
