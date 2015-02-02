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
using namespace std;

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
   string _attr_val;
   GffAttr()=default;
   GffAttr(int id, string val){
      _attr_id = id;
      _attr_val = val;
   }
};

class GffInfo{
   //Info is a id->name type
public:
   int _id;
   string _name;
   GffInfo() = default;
   GffInfo(int id, string name){
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
   vector<GffInfo> _gff_info_vec;
public:
   unordered_map<string, int> _name2id;
   int addInfo(const string name);
};

class GffInfoTable{
public:
   unique_ptr<GffInfoVec> _seq_names;
   unique_ptr<GffInfoVec> _attr_names;
   GffInfoTable(){
      _seq_names = unique_ptr<GffInfoVec>(new GffInfoVec());
      _attr_names = unique_ptr<GffInfoVec>(new GffInfoVec());
   }
};

class GffLine{
   char* _info;
   void extractAttr(const string attr, string &val);
   bool _is_gff3;
   char* _line;
public:
   char* _dupline;
   int _llen;
   string _chrom;
   string _source;
   char _strand;
   string _gffline_type; // the type indicated by gff line
   uint _start;
   uint _end;
   float _score;
   bool _skip;
   //bool _is_cds; // for future
   GffFeat_t _feat_type; // the type which has been parsed
   char _phase;
   string _ID;
   string _name;
   string _parent;
   vector<string> _parents;
   GffLine(const char* l);
   GffLine(const GffLine &l);
   GffLine(GffLine &&l);
   GffLine& operator = (const GffLine &rhs) = delete;
   GffLine& operator = (GffLine &&rhs) = default;
   GffLine();
   ~GffLine();
};
typedef shared_ptr<GffLine> LinePtr;

class GffReader;
class GffObj{
protected:
   GenomicInterval _iv;
   int _seq_id;
   vector<GffAttr> _attrs; // Attr id is the vector idx
   double _score;
   char _phase;
   string _source;
   GffReader & _greader;
public:
   static unique_ptr<GffInfoTable> _infotable;
   GffObj() = default;
   GffObj(LinePtr gl, GffReader & greader);
   virtual ~GffObj(){};
   GffObj(const GffObj &other) = default;
   GffObj(GffObj &&other) = default;
   GffObj& operator=(const GffObj &rhs) = default;
   GffObj& operator=(GffObj &&rhs) = default;
   virtual int seq_id() const { return _seq_id;}
   virtual const char strand() const { return _iv.strand();}
};
class GffmRNA;

class GffExon: public GffObj{
   vector<string> _parent_mrnas;
   string _parent_gene;
public:
   bool _within_3UTR;
   bool _within_5UTR;
   vector<GffmRNA*> _parents;
   string _exon_id;
   string _exon_name;
   GffExon(LinePtr gl, GffReader & greader);
   const string parent_mrna() const { return _parent_mrna;}
};

class GffIntron: public GffObj{
   vector<string> _parent_mrnas;
   string _parent_gene;
public:
   bool _within_3UTR;
   bool _within_5UTR;
   int _intron_id;
   string _intron_name;
   vector<GffmRNA*> _parents;
   GffIntron(LinePtr gl, GffReader & greader);
   const string parent_mrna() const { return _parent_mrna;}
};

class GffUTR: public GffObj{
   string _parent_mrna;
public:
   int _utr_id;
   string _utr_name;
   GffmRNA* _parent;
   GffUTR(LinePtr gl, GffReader & greader);
   const string parent_mrna() const { return _parent_mrna;}
};

class GffLoci: public GffObj{
public:
   string _gene_id;
   string _gene_name;
   vector<GffmRNA* > _mrnas;
   vector<GffExon> _non_dup_exons;
   vector<GffIntron> _non_dup_introns;
   GffLoci(LinePtr gl, GffReader & greader);
};

class GffmRNA: public GffObj{
   string _parent_gene;
public:
   string _transcript_id;
   string _transcript_name;
   GffLoci* _parent;
   vector<*GffExon > _exons;
   vector<*GffIntron > _introns;
   vector<GffUTR> _UTRs;
   GffmRNA(LinePtr gl, GffReader & greader);
   const string parent_gene() const { return _parent_gene;}
};


class GffSeqData {
public:
   vector<GffmRNA> _forward_rnas;
   vector<GffmRNA> _reverse_rnas;
   vector<GffmRNA> _unstranded_rnas;
   vector<GffLoci> _genes;
   GffmRNA& last_f_rna()  { return _forward_rnas.back();}
   GffmRNA& last_r_rna()  { return _reverse_rnas.back();}
   GffmRNA& last_u_rna() { return _unstranded_rnas.back();}
   GffLoci& last_gene()  { return _genes.back();}
};


class GffReader: public SlineReader{
   string _fname;
public:
   vector<unique_ptr<GffSeqData> > & _g_seqs;
   LinePtr _gfline;
   GffReader(vector<unique_ptr<GffSeqData>> & gseqs,const char* f=NULL);
   bool nextGffLine();
   void readAll();
};

#endif

