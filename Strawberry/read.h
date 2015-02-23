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

   CigarOpCode opcode : 3;
   uint32_t length : 29;

   CigarOp(CigarOpCode o, uint32_t l) : opcode(o), length(l) {}
   bool operator==(const CigarOp& rhs) const { return opcode == rhs.opcode && length == rhs.length; }

};

typedef uint64_t ReadID;
typedef int RefID;
class ReadHit{
private:

   ReadID _read_id;
   GenomicInterval _iv;
   RefID _partner_ref_id;
   int _partner_pos;
   int _num_mismatch = -1;
   int _num_hit = 1;
   int _trans_left=0; // position relate to transcriptom
   uint32_t _sam_flag = 0; //bitwise FLAG
   float _read_mass = 0.0;

public:
   vector<CigarOp> _cigar;
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
   int partner_pos() const;
   RefID ref_id() const; // chromosome id or scaffold id containing the read
   int num_mismatch() const;
   bool is_singleton() const;
   uint left() const;
   uint right() const;
   char strand() const;
   double mass() const;
   //vector<CigarOp> cigars() const;
};


// Now it is used only to convert read name to read id
//
class ReadTable
{
public:
   // This function should NEVER return zero
   ReadID get_id(const string& name);

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
protected:
   static const unsigned MAX_HEADER_LEN = 4 * 1024 * 1024; // 4 MB
   static const size_t kHitBufMaxSize = 10 * 1024;
   ReadTable& _reads_table;
   char _hit_buf[kHitBufMaxSize];
   int _num_seq_header_recs = 0;
public:
   RefSeqTable& _ref_table;
   HitFactory(ReadTable &reads_table, RefSeqTable &ref_table);
   HitFactory(HitFactory &rhs) = delete; //non-copible class.
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

//PairedHit solely own the ReadHit objects by
//using unique_ptr
class PairedHit{
   ReadHitPtr _left_read ;
   ReadHitPtr _right_read;
   double _collapse_mass = 0.0;

public:
   PairedHit() = default;
   PairedHit(ReadHitPtr leftRead, ReadHitPtr rightRead);
   const ReadHitPtr left_read() const;
   void set_left_read(ReadHitPtr lr);
   char strand() const;
   const ReadHitPtr right_read() const;
   void set_right_read(ReadHitPtr rr);
   bool is_paired() const;
   RefID ref_id() const;
   uint left_pos() const;
   uint right_pos() const;
   uint edit_dist() const;
   bool contains_splice() const;
   ReadID read_id() const;

   bool paried_hit_lt(const PairedHit &rhs) const;
};


#endif /* STRAWB_READ_H_ */
