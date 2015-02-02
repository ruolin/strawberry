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
   int addInfo(string name);
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
typedef unique_ptr<GffLine> LinePtr;

class GffObj{
protected:
   GenomicInterval _iv;
   static unique_ptr<GffInfoTable> _infotable;
   int _seq_id;
   vector<GffAttr> _attrs;
   double _score;
   char _phase;
   char _souce;
public:
   GffObj() = default;
   GffObj(LinePtr gl);
   virtual ~GffObj(){};
   GffObj(const GffObj &other) = default;
   GffObj(GffObj &&other) = default;
   GffObj& operator=(const GffObj &rhs) = default;
   GffObj& operator=(GffObj &&rhs) = default;
};
class GffmRNA;

class GffExon: public GffObj{
public:
   bool _within_3UTR;
   bool _within_5UTR;
   shared_ptr<GffmRNA> _parent;
   int _exon_id;
   string _exon_name;
   GffExon(LinePtr gl);
};

class GffIntron: public GffObj{
public:
   bool _within_3UTR;
   bool _within_5UTR;
   int _intron_id;
   string _intron_name;
   shared_ptr<GffmRNA> _parent;
   GffIntron(LinePtr gl);
};

class GffUTR: public GffObj{
public:
   int _utr_id;
   string _utr_name;
   shared_ptr<GffmRNA> _parent;
   GffUTR(LinePtr gl);
};

class GffLoci: public GffObj{
public:
   int _gene_id;
   string _gene_name;
   vector<shared_ptr<GffmRNA> > _mrnas;
   set<GffExon> _non_dup_exons;
   set<GffIntron> _non_dup_introns;
   GffLoci(LinePtr gl);
};

class GffmRNA: public GffObj{
public:
   int _transcript_id;
   string _transcript_name;
   shared_ptr<GffLoci> _parent;
   vector<shared_ptr<GffExon> > _exons;
   vector<shared_ptr<GffIntron> > _introns;
   vector<GffUTR> _UTRs;
   GffmRNA(LinePtr gl);
};





class GffReader: public SlineReader{
   string _fname;
public:
   LinePtr _gfline;
   GffReader(const char* f=NULL);
   bool nextGffLine();
   void readAll();
};

#endif

