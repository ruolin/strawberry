/*
 * read_start.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: Ruolin Liu
 */

#ifndef GENOME_H
#define GENOME_H

#include <vector>
#include <unordered_map>
#include <memory>
#include <set>
#include "contig.h"
#include "common.h"

enum GffFeat_t{
   OTHERS = 0,
   GENE,
   mRNA,
   EXON,
   INTRON,
   UTR,
   CDS,
   STOP_CODON,
   START_CODON
};

class GffAttr{
public:
   int _attr_id;
   std::string _attr_val;
   GffAttr()=default;
   GffAttr(int id, std::string val){
      _attr_id = id;
      _attr_val = val;
   }
};

class GffInfo{
   /*
    * Info is a id->name pair.
    * It is intended to be use in more general case.
    * But right now only use for ref_id -> ref_seq_name
    */
public:
   int _id;
   std::string _name;
   GffInfo() = default;
   GffInfo(int id, std::string name){
      _id = id;
      _name = name;
   }
   bool operator==(GffInfo &rhs){
      return (_name.compare(rhs._name) == 0);
   }
   bool operator<(GffInfo &rhs){
      return (_name.compare(rhs._name) < 0 );
   }
};

class GffInfoVec{
   std::vector<GffInfo> _gff_info_vec;
public:
   std::unordered_map<std::string, int> _name2id;
   int addInfo(const std::string name);
};

class GffInfoTable{
   /*
    * Lookup table for reference sequence names and attribute names
    * FIXME: Attribute names has not been used yet.
    */
public:
   std::unique_ptr<GffInfoVec> _seq_names;
   //std::unique_ptr<GffInfoVec> _attr_names;
   GffInfoTable(){
      _seq_names = std::unique_ptr<GffInfoVec>(new GffInfoVec());
     // _attr_names = std::unique_ptr<GffInfoVec>(new GffInfoVec());
   }
};

class GffLine{
   char* _info;
   void extractAttr(const std::string attr, std::string &val);
   bool _is_gff3;
   char* _line;
public:
   char* _dupline;
   int _llen;
   uint _end;
   uint _start;
   std::string _source;
   Strand_t _strand;
   bool _skip;
   float _score;
   std::string _chrom;
   std::string _gffline_type; // the type indicated by gff line
   //bool _is_cds; // for future
   GffFeat_t _feat_type; // the type which has been parsed
   char _phase;
   std::string _ID;
   std::string _name;
   std::string _parent;
   std::vector<std::string> _parents;
   GffLine(const char* l);
   GffLine(const GffLine &l);
   GffLine(GffLine &&l);
   GffLine& operator = (const GffLine &rhs) = delete;
   GffLine& operator = (GffLine &&rhs) = default;
   const std::string parent() const
   {
      if(_parents.size() != 1){
         std::cerr<<"No parent gene or multiple parents for a mRNA object in "<< _dupline<<std::endl;;
         exit(0);
      }

      return _parents[0];
   }
   GffLine();
   ~GffLine();
};
typedef std::shared_ptr<GffLine> LinePtr;

class GffReader;
class GffObj{
protected:
   double _score;
   std::vector<GffAttr> _attrs; // Attr id is the vector idx
   char _phase;
   std::string _source;
   GffReader & _greader;
public:
   GenomicInterval _iv;
   static std::unique_ptr<GffInfoTable> _infotable;
   GffObj() = default;
   GffObj(LinePtr gl, GffReader & greader);
   virtual ~GffObj(){};
   GffObj(const GffObj &other) = default;
   GffObj(GffObj &&other) = default;
   GffObj& operator=(const GffObj &rhs) = default;
   GffObj& operator=(GffObj &&rhs) = default;
   virtual int seq_id() const { return _iv.seq_id();}
   virtual Strand_t strand() const { return _iv.strand();}
   virtual bool operator==(const GffObj &rhs){
      return _iv == rhs._iv;
   }
   virtual bool operator<(const GffObj &rhs){
      if(_iv.strand() == rhs._iv.strand()){
         if(_iv.strand() == Strand_t::StrandMinus) return _iv > rhs._iv;
         else return _iv < rhs._iv;
      }
      return false;
   }
};

class GffmRNA;
class GffLoci;



class GffExon: public GffObj{
   GffLoci* const _parent_gene;
public:
   bool _within_3UTR;
   bool _within_5UTR;
   std::vector<GffmRNA*> _parent_mrnas;
   std::string _exon_id;
   std::string _exon_name;
   GffExon(LinePtr gl, GffmRNA* mrna, GffLoci* const gene, GffReader & greader);
   GffLoci* const parent_gene()
   {
      return _parent_gene;
   }
   void add_parent_mrna(const std::vector<GffmRNA*> & mrnas){
      for(size_t i=0; i<mrnas.size(); i++){
         _parent_mrnas.push_back(mrnas[i]);
      }
   }
};



/*
 * I decided not to use GffIntron for now but reserved for the future use.
 */
class GffIntron: public GffObj{
   const GffLoci* const _parent_gene;
public:
   bool _within_3UTR;
   bool _within_5UTR;
   std::string _intron_id;
   std::string _intron_name;
   std::vector<GffmRNA*> _parents_mrnas;
   GffIntron(LinePtr gl, GffLoci &gene, GffReader & greader);
   const GffLoci* const parent_gene() const
   {
      return _parent_gene;
   }
};



class GffUTR: public GffObj{
   const GffLoci* const _parent_gene;
public:
   std::string _utr_id;
   std::string _utr_name;
   std::vector<GffmRNA*> _parents_mrnas;
   GffUTR(LinePtr gl, GffReader & greader);
   const GffLoci* const parent_gene() const
   {
      return _parent_gene;
   }
};

class GffCDS: public GffObj{
   const GffLoci* const _parent_gene;
public:
   std::string _cds_id;
   std::string _cds_name;
   std::vector<GffmRNA*> _parents_mrnas;
   GffCDS(LinePtr gl, GffmRNA mrna, GffLoci &gene, GffReader &greader);
   const GffLoci* const parent_gene() const
   {
      return _parent_gene;
   }
};

using exonPtr = std::unique_ptr<GffExon>;
using intronPtr = std::unique_ptr<GffIntron>;
using utrPtr= std::unique_ptr<GffUTR>;
using cdsPtr = std::unique_ptr<GffCDS>;

class GffLoci: public GffObj{
public:
   std::vector<GffmRNA* > _mrnas;
   std::string _gene_id;
   std::string _gene_name;
   std::vector<exonPtr> _non_dup_exons;
   std::vector<intronPtr> _non_dup_introns;
   std::vector<utrPtr> _non_dup_utrs;
   GffLoci(LinePtr gl, GffReader & greader);
   void add_mRNA(GffmRNA *gffmrna)
   {
      _mrnas.push_back(gffmrna);
   }
   void add_exon(exonPtr exon, GffmRNA* exon_parent);
   GffmRNA* getRNA(const std::string rna);
   int num_mRNAs() const{
      return _mrnas.size();
   }
};


class GffmRNA: public GffObj{
   GffLoci* const _parent;
public:
   std::string _transcript_id;
   std::string _transcript_name;
   std::vector<GffExon* > _exons;
   std::vector<GffIntron* > _introns;
   std::vector<GffUTR* > _UTRs;
   GffmRNA(LinePtr gl, GffLoci* gene, GffReader & greader);
   GffLoci* const getParentGene() const{
      return _parent;
   }
   void add_exon(GffExon* exon){
      _exons.push_back(exon);
   }
};

using mrnaPtr = std::unique_ptr<GffmRNA>;
using genePtr = std::unique_ptr<GffLoci>;

class GffTree {
    /*
     *
     */
public:
   std::vector<mrnaPtr> _forward_rnas;
   std::vector<mrnaPtr> _reverse_rnas; // in each GffmRNA object, exon order is from small-to-large
   std::vector<mrnaPtr> _unstranded_rnas;
   std::vector<genePtr> _genes; // in each GffLoci object, non_dup_exons contain the exons
                           // according to the order in gff file. Usually in Gff3 format
                           // the order for minus strand is large-to-small.
   std::string _g_seq_name;

   explicit GffTree(const std::string &g_seq_name):
         _g_seq_name(g_seq_name)
   {}
   GffTree() = default;

   int get_gseq_id() const{
      assert(!_genes.empty());
      return _genes.front()->_iv.seq_id();
   }

   void set_gseq_id(int id)
   {
      for(auto &g : _genes){
         g->_iv.set_seq_id(id);
         for(auto &e : g->_non_dup_exons){
            e->_iv.set_seq_id(id);
         }
      }
      for(auto &m : _forward_rnas){
         m->_iv.set_seq_id(id);
      }
      for(auto &m: _reverse_rnas){
         m->_iv.set_seq_id(id);
      }
      for(auto &m: _unstranded_rnas){
         m->_iv.set_seq_id(id);
      }
   }

   GffmRNA* last_f_rna()
   {
      return &(*_forward_rnas.back());
   }

   GffmRNA* last_r_rna()
   {
      return &(*_reverse_rnas.back());
   }

   GffmRNA* last_u_rna()
   {
      return &(*_unstranded_rnas.back());
   }

   GffLoci* last_gene()
   {
      return &(*_genes.back());
   }

   void addGene(genePtr gene)
   {
      _genes.push_back(move(gene));
   }

   void addPlusRNA(mrnaPtr mrna)
   {
      _forward_rnas.push_back(move(mrna));
   }

   void addMinusRNA(mrnaPtr mrna)
   {
      _reverse_rnas.push_back(move(mrna));
   }

   void addUnstrandedRNA(mrnaPtr mrna)
   {
      _unstranded_rnas.push_back(move(mrna));
   }

   GffLoci* findGene(const std::string gene_id);
   GffmRNA* findmRNA(const std::string mrna_id, const Strand_t strand);

   std::vector<Contig> loadRnas() const;
};


class GffReader: public SlineReader{
   std::string _fname;
public:
   std::vector<std::unique_ptr<GffTree> >  _g_seqs;
   LinePtr _gfline;
   GffReader(const char* fname, FILE* stream=NULL);
   void load_gff(const char* f=NULL);
   bool nextGffLine();
   void readAll();
   void addGseq(std::unique_ptr<GffTree> gseq){
      _g_seqs.push_back(move(gseq));
   }
   void reverseExonOrderInMinusStrand();
};

#endif



