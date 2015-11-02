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
#include<iostream>
#include "common.h"

class PairedHit;
class ReadHit;
class RefSeqTable;
class ExonBin;
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
   double _avg_cov; // average coverage
   GenomicFeature(const Match_t& op, uint offset, int len);
   void left(uint left);
   uint left() const;
   void right(uint right);
   uint right() const;
   void avg_doc(double coverage);
   double avg_doc() const;

   static bool overlaps(const GenomicFeature& lhs, const GenomicFeature& rhs);
   static int overlap_len(const GenomicFeature &lhs, const GenomicFeature &rhs);
   static bool overlap_in_genome(const GenomicFeature& feat, uint s, uint e);

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
public:
   std::vector<GenomicFeature> _genomic_feats;
   bool _is_ref;
   //Contig() = default;
   Contig(const PairedHit& ph);
   Contig(RefID ref_id, Strand_t strand, const std::vector<GenomicFeature> &feats, double mass,bool is_ref);
   Contig(const ExonBin& eb);
   const std::string annotated_trans_id() const;
   void annotated_trans_id(std::string str);
   uint left() const;
   uint right() const;
   static bool overlaps_only_on_exons(const Contig &ct, const GenomicFeature & gf);
   static bool overlaps_directional(const Contig &lhs, const Contig &rhs);
   static bool is_contained_in(const Contig &small, const Contig &large);
   static int infer_frag_len(const Contig &isoform, const Contig &hit);
   bool operator<(const Contig &rhs) const;
   RefID ref_id() const;
   Strand_t strand() const;
   size_t featSize() const;
   float mass() const;
   void mass(float m);
   void print2gtf(FILE *pFile, const RefSeqTable &ref_lookup, int gene_id, int tscp_id);
};

class Isoform{

public:
   Contig _contig;
   string _gene_id;
   string _isoform_id;
   double _bais_factor;
   //Isoform() = default;
   Isoform(Contig contig, string gene, string isoform);
};

class ExonBin{
   /*
    * ExonBin is a set a continuous exons defined by an overlapping fragment.
    */
private:
   //uint _left_most;
   //uint _right_most;
   int _read_num;
   RefID _ref_id;
public:
   std::set<uint> _coords;
   std::vector<const GenomicFeature*> _exon_in_bin;
   std::vector<string> _parent_isoforms;
   std::vector<const Contig*> _hits;
   ExonBin(RefID ref_id);
   ExonBin(const ExonBin &rhs) = delete;
   ExonBin& operator=(const ExonBin &rhs) = delete;
   ExonBin(ExonBin &&rhs) = default;
   ExonBin& operator=(ExonBin &&rhs) = default;
   uint left_most() const;
   uint right_most() const;
   bool insert_exon(const GenomicFeature *exon);
   void read_num_increase_by_1();
   int read_num() const;
   void add_hit(const Contig* hit);
   void add_hit(const ExonBin& exb);
   void add_isoform(const vector<Isoform> &isoforms);
   RefID ref_id() const;
};

#endif /* CONTIG_H_ */
