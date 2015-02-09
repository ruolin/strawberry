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
#include "read.h"

enum Op_t {S_MATCH, S_INTRON, S_UNKNOWN};

struct GenomicFeature{
   Op_t _code;
   int _genomic_offset;
   int _genomic_length;
   GenomicFeature(const Op_t& cc, int offset, int len);

   void g_left(int left);
   int g_left() const;
   void g_right(int right);
   int g_right() const;

   static bool overlap_in_genome(const GenomicFeature& lhs, const GenomicFeature& rhs);

   bool contains(const GenomicFeature& other) const;

   bool properly_contains(const GenomicFeature& other) const;

   static int match_length(const GenomicFeature &op, int left, int right);

   bool operator==(const GenomicFeature & rhs) const;

   bool operator<(const GenomicFeature & rhs) const;
};



/* Reduced representation of GffmRNA using a vector of GenomicFeatures which include intron
 * exon, and unknown*/

class Contig{
   std::vector<GenomicFeature> _genomic_feats;
   bool _is_ref;
   std::vector<const PairedHit* > _mates_in_contig;
   std::string _seq;
   std::string _annotated_gene_id;
   std::string _annotated_trans_id;
   RefID _ref_id;
   char _strand;
public:
   Contig(RefID ref_id, char strand, vector<GenomicFeature> &ops, bool is_ref);
};

#endif /* CONTIG_H_ */
