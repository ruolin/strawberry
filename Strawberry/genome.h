/*
 * read_start.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: Ruolin Liu
 */

#ifndef GENOME_H
#define GENOME_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <memory>
#include "GBase.h"
static const char strand_plus = '+';
static const char strand_minus = '-';
static const char strand_unknown = '.';

#define GFF_LINELEN 2048
// to do:
enum gffExonType { exGffUnknown = 0, exGffUTR = 1, exGffStop = 2,
   exGffStart =3, exGffCDS = 4, exGffExon=5};

using namespace std;

class GenomicInterval {
   uint _start; //start<end always!
   uint _end;
   char* _chrom;
   char _strand;
 public:
   GenomicInterval();
   GenomicInterval(const char* chr,uint s = 0,uint e = 0, char o = '.');
   GenomicInterval( const GenomicInterval& other);
   ~GenomicInterval(){
      //free(chrom);
      GFREE(_chrom);
   }

  uint getStart();
  uint getEnd();
  char getStrand();
  char* getChrom();


  //check for overlap with other segment
  uint len();
  bool overlap(GenomicInterval* d, bool nonStrandness = true);

  bool isContainedIn(GenomicInterval *d, bool nonStrandness = true);

  bool contain(GenomicInterval *d, bool nonStrandness = true);

  //return the length of overlap between two segments
  uint overlapLen(GenomicInterval* r) ;
  //comparison operators required for sorting
  bool operator==(GenomicInterval& d){
     if ( strcmp(d._chrom, _chrom)) return false;
     if( d._strand != strand_unknown && _strand != strand_unknown && d._strand != _strand) return false;
     return (_start==d._start && _end==d._end);
  }
  bool operator>(GenomicInterval& d){
     if ( strcmp(d._chrom, _chrom)) GError("cannot compare for different chrom");
     return (_start==d._start)?(_end>d._end):(_start>d._start);
  }
  bool operator<(GenomicInterval& d){
     if ( strcmp(d._chrom, _chrom)) GError("cannot compare for different chrom");
     return (_start==d._start)?(_end<d._end):(_start<d._start);
  }
};


class GffAttr{
   public:
   int _attr_id;
   char* _attr_val;
   GffAttr(int an_id, const char* av=NULL);
  ~GffAttr();
  void setAttrValue(const char *av);
};

class GffExon{
public:
   GenomicInterval _exon_iv;
   string _exon_id={}; // optional field if gff file has such information
   int _exon_num; // the exon apparence in gff file.
   string _parent_transcript_id;
   string _parent_gene_id={};
   double _score; // gff score column
   char _phase;
   int _tstart=0; // relative starts in mRNA
   int _tend=0; // relative ends in mRNA
   GffExon( GenomicInterval iv, int exnum, string tid,
         double score, char phase);
   ~GffExon();
};

typedef shared_ptr<GffExon> exonPtr;

class GffTranscript{
public:
   vector<exonPtr> _exons; // transcripts
   GenomicInterval _trans_iv;
   string _parent_gene_id;
   string _trans_id;
   GffTranscript() = default;
   GffTranscript(GenomicInterval iv, string tid);
   ~GffTranscript(){_exons.clear();}
   bool setParentGeneID(string gid);
   void addExon(shared_ptr<GffExon> e);
};
typedef shared_ptr<GffTranscript> transPtr;

class GffGene{
public:
   unordered_map<string, transPtr> _transcripts;
   unordered_map<int, exonPtr> _unique_exons;
   GenomicInterval _gene_iv;
   string _gene_id;
   unsigned short _transcript_num;
   GffGene(GenomicInterval iv, string _gene_id);
   ~GffGene(){ _transcripts.clear(); _unique_exons.clear();}
   void addTranscript(shared_ptr<GffTranscript> &t);

};
typedef shared_ptr<GffGene> geneSharedPtr;
typedef unique_ptr<GffGene> geneUniqPtr;

class GffLine{
public:
   char* _dupline;
   char* _line;
   int _llen;
   char* _chrom;
   char* _track;
   char _strand;
   char* _ftype;
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
   char* _ID;
   char* _name;
   char** _parents;
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

