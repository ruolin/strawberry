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

#define GFF_LINELEN 2048
using namespace std;

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


class GffObj{
protected:
   GenomicInterval _iv;
   static unique_ptr<GffInfoTable> _infotable = nullptr;
   int _seq_id;
   vector<GffAttr> _attrs;
   double _score;
   char _phase;
   char _souce;
public:
   GffObj() = default;
   GffObj(GffLine gl);
   virtual ~GffObj(){};
   GffObj(const GffObj &other) = default;
   GffObj(GffObj &&other) = default;
   GffObj& operator=(const GffObj &rhs) = default;
   GffObj& operator=(GffObj &&rhs) = default;
};

class GffExon: public GffObj{
public:
   bool _within_3UTR;
   bool _within_5UTR;
   shared_ptr<GffmRNA> _parent;
   int _exon_id;
   string _exon_name;
   GffExon(GffLine gl);
};

class GffmRNA: public GffObj{
public:
   int _transcript_id;
   string _transcript_name;
   shared_ptr<GffLoci> _parent;
   vector<shared_ptr<GffExon> > _exons;
   vector<shared_ptr<GffIntron> > _introns;
   vector<Gff5UTR> _5UTRs;
   vector<Gff3UTR> _3UTRs;
   GffmRNA(GffLine gl);
};

class GffLoci: public GffObj{
public:
   int _gene_id;
   string _gene_name;
   vector<shared_ptr<GffmRNA> > _mrnas;
   set<GffExon> _non_dup_exons;
   set<GffIntron> _non_dup_introns;
   GffLoci(GffLine gl);
};

class GffIntron: public GffObj{
public:
   bool _within_3UTR;
   bool _within_5UTR;
   int _intron_id;
   string _intron_name;
   shared_ptr<GffmRNA> _parent;
   GffIntron(GffLine gl);
};

class Gff5UTR: public GffObj{
public:
   int _5utr_id;
   string _5utr_name;
   shared_ptr<GffmRNA> _parent;
   Gff5UTR(GffLine gl);
};

class Gff3UTR: public GffObj{
public:
   int _3utr_id;
   string _3utr_name;
   shared_ptr<GffmRNA> _parent;
   Gff3UTR(GffLine gl);
};

class GffLine{
public:
   char* _dupline;
   char* _line;
   int _llen;
   string _chrom;
   string _source;
   char _strand;
   string _ftype;
   char* _info;
   uint _start;
   uint _end;
   double _score;
   bool _skip;
   bool _is_gff3;
   //bool _is_cds; // for future
   bool _is_exon;
   int _exontype;
   bool _is_transcript;
   bool _is_gene;
   char _phase;
   string _ID;
   string _name;
   vector<string> _parents;
   int _parents_len;
   int _num_parents;
   GffLine(const char* l);
   void discardParent();
   GffLine(GffLine *l);
   GffLine();
   char* extractAttr(const char* pre,
         bool caseStrict=false,
         bool enforce_GTF2=false);
   ~GffLine();
};
typedef shared_ptr<GffLine> LinePtr;


class GffReader{
   friend class GffTranscript;
   friend class GffExon;
   friend class GffLine;
   friend class GffGene;
   char* _fname;
   char* _linebuf;
   int _buflen;
   off_t _fpos;
protected:
   FILE *_fh;
public:
   unordered_map<string, geneSharedPtr> _genes={};
   unordered_map<string, string> _trans2gene={};
   LinePtr _gffline;
   GffReader(char* f=NULL);
   ~GffReader();
   bool nextGffLine();
   void readAll();
};

#endif

