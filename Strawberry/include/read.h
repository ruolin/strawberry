/*
 * read.h
 *
 *  Created on: Nov 3, 2014
 *      Author: RUOLIN
 */

#ifndef STRAWB_READ_H_
#define STRAWB_READ_H_
#include<string>
#include<bam/sam.h>
#include<vector>
#include<memory>
#include<string>
#include<unordered_map>
#include<cassert>
#include "common.h"
#ifdef DEBUG
   #include <iostream>
#endif
using namespace std;
/*
 * COMMON PARAMETERS
 */

enum CigarOpCode
   {
      MATCH = BAM_CMATCH,
      INS = BAM_CINS,
      DEL = BAM_CDEL,
      REF_SKIP = BAM_CREF_SKIP,
      SOFT_CLIP = BAM_CSOFT_CLIP,
      HARD_CLIP = BAM_CHARD_CLIP,
      PAD = BAM_CPAD,
      MISMATCH = 7
};

struct CigarOp
{

   CigarOpCode _type : 3;
   uint32_t _length : 29;

   CigarOp(CigarOpCode o, uint32_t l) : _type(o), _length(l) {}
   bool operator==(const CigarOp& rhs) const { return _type == rhs._type && _length == rhs._length; }

};

class ReadHit{
private:

   ReadID _read_id;
   GenomicInterval _iv;
   vector<CigarOp> _cigar;
   RefID _partner_ref_id;
   uint _partner_pos;
   int _num_mismatch = -1;
   int _num_hits = 1;
   int _trans_left=0; // position relate to transcriptom
   uint32_t _sam_flag = 0; //bitwise FLAG
   float _read_mass = 0.0;
public:
   ReadHit() = default;
   ReadHit( ReadID readID,
         GenomicInterval iv,
         const vector<CigarOp> & cigar,
         RefID partnerRef,
         int partnerPos,
         int numMismatch,
         int numHit,
         uint32_t samFlag,
         float mass
         );

   uint read_len() const;
   ReadID read_id() const;
   bool contains_splice() const;
   GenomicInterval interval() const;
   RefID partner_ref_id() const;
   const vector<CigarOp>& cigar() const;
   uint partner_pos() const;
   RefID ref_id() const; // chromosome id or scaffold id containing the read
   int num_mismatch() const;
   int numHits() const;
   bool is_singleton() const;
   uint left() const;
   uint right() const;
   Strand_t strand() const;
   float mass() const;
   bool is_first() const;
   bool reverseCompl() const;
   bool operator==(const ReadHit& rhs) const; // not considering read orientation
   bool operator!=(const ReadHit& rhs) const;
   //vector<CigarOp> cigars() const;
};


// Now it is used only to convert read name to read id
//
class ReadTable
{

public:
   // This function should NEVER return zero
   ReadID get_id(const string& name);
   uint _read_len_abs = 0;
   vector<int> _frag_dist;
private:
   // This is FNV-1, see http://en.wikipedia.org/wiki/Fowler_Noll_Vo_hash
   inline uint64_t hashString(const char* __s)
   {
      uint64_t hash = 0xcbf29ce484222325ull;
      for ( ; *__s; ++__s)
      {
         hash *= 1099511628211ull;
         hash ^= *__s;
      }
      return hash;
   }
};

class InsertSize{
public:
   double _total_reads;
   vector<double> _emp_dist;
   int _start_offset;
   int _end_offset;
   double _mean;
   double _sd;
   InsertSize();
   InsertSize(const vector<int> frag_lens);
   double emp_dist_pdf(uint insert_size) const;
   //double truncated_normal_pdf(uint insert_size) const;
};


class RefSeqTable
{

private:
//_id is the observation order.
//_id start from 0 which is used as the index in vector<GffSeqData>;
   string _seq;
   bool _keep_seq;
   unordered_map<string, int> _name2id;
   vector<string> _id2name;

public:
   RefSeqTable(bool keep_seq) : _keep_seq(keep_seq){}
   int get_id(const string& name);
   int size() const{
      return _name2id.size();
   }
   const string ref_name(int id) const{
      return _id2name[id];
   }
};


class HitFactory
{
private:
   virtual platform_t str2platform(const string pl_str);
protected:
   static const unsigned MAX_HEADER_LEN = 4 * 1024 * 1024; // 4 MB
   static const size_t kHitBufMaxSize = 10 * 1024;
   char _hit_buf[kHitBufMaxSize];
   int _num_seq_header_recs = 0;
   AssayProperties _assay_props;
public:
   ReadTable& _reads_table;
   RefSeqTable& _ref_table;
   HitFactory(ReadTable &reads_table, RefSeqTable &ref_table);
   HitFactory(const HitFactory &rhs) = delete; //non-copible class.
   HitFactory& operator=(const HitFactory &rhs) = delete; // non-copible class
   HitFactory(HitFactory &&rhs) = default;
   HitFactory& operator=(HitFactory &&rhs) = default;
   virtual ~HitFactory() = default;
   virtual bool recordsRemain() const = 0;

   // return false when no record or EOF found
   virtual bool nextRecord(const char* &buf, size_t& buf_size) = 0;
   virtual bool getHitFromBuf(const char* bwt_buf, ReadHit& bh)=0;
   virtual RefSeqTable& ref_table() { return _ref_table; }
   virtual ReadTable& reads_table(){return _reads_table;}
   virtual void undo_hit() = 0;
   virtual bool parse_header_line(const string& hline);
   virtual bool inspect_header() = 0;
   virtual void reset() = 0;
   virtual const AssayProperties& assay_properties();
};

class BAMHitFactory : public HitFactory
{
private:

   samfile_t* _hit_file;
   int64_t _curr_pos;
   int64_t _beginning;

   bam1_t _next_hit;
   bool _eof_encountered;

public:
   BAMHitFactory(const string& bam_file_name,
                 ReadTable& read_table,
                 RefSeqTable &ref_table) throw(runtime_error);
   ~BAMHitFactory();
   bool recordsRemain() const;
   void markCurrPos();
   void reset();
   void undo_hit();
   bool nextRecord(const char* &buf, size_t& buf_size);
   bool getHitFromBuf(const char* bwt_buf, ReadHit& bh);
   bool inspect_header();
};

typedef shared_ptr<ReadHit> ReadHitPtr;


class PairedHit{
   double _collapse_mass = 0.0;
public:
   ReadHitPtr _right_read ;
   ReadHitPtr _left_read ;
   PairedHit() = default;
   PairedHit(ReadHitPtr leftRead, ReadHitPtr rightRead);
   PairedHit& operator=(const PairedHit& rhs){
      _left_read = rhs._left_read;
      _right_read = rhs._right_read;
      _collapse_mass = rhs._collapse_mass;
   }
   const ReadHit& left_read_obj() const;
   void set_left_read(ReadHitPtr lr);
   Strand_t strand() const;
   const ReadHit& right_read_obj() const;
   void set_right_read(ReadHitPtr rr);
   bool is_paired() const;
   RefID ref_id() const;
   uint left_pos() const;
   uint right_pos() const;
   uint edit_dist() const;
   bool contains_splice() const;
   ReadID read_id() const;
   int numHits() const;
   bool is_multi() const;
   float raw_mass() const;
   bool operator==(const PairedHit& rhs) const;
   bool operator!=(const PairedHit& rhs) const;
   bool operator<(const PairedHit& rhs) const;
   bool paried_hit_lt(const PairedHit &rhs) const;
   void add_2_collapse_mass(float add);
   float collapse_mass() const;
};

void mean_and_sd_insert_size(const vector<int> & vec, double & mean, double &sd);
#endif /* STRAWB_READ_H_ */
