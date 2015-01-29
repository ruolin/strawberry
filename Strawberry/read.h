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
#include<unordered_map>
#include "contig.h"
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

class ReadHit{
private:
   static const int max_partner_dist = 50000;

   ReadID _read_id;
   GenomicInterval _iv;
   string _partner_ref;
   int _partner_pos;
   vector<CigarOp> _cigar;
   int _num_mismatch = -1;
   int _num_hit = 1;
   int _trans_left=0; // position relate to transcriptom
   uint32_t _sam_flag = 0; //bitwise FLAG
   float _read_mass = 0.0;

public:
   ReadHit() = default;
   ReadHit( ReadID readID,
         GenomicInterval iv,
         const vector<CigarOp> & cigar,
         string partnerRef,
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
   string partner_ref() const;
   int partner_pos() const;
   string ref() const; // chromosome or scaffold containing the read
   int num_mismatch() const;
   bool is_singleton() const;
   int left() const;
   int right() const;
   char strand() const;
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

   struct SequenceInfo{
      uint32_t _observation_order;
      unique_ptr<DNABitSet> _seq;
      SequenceInfo(uint32_t order,
                   unique_ptr<DNABitSet> seq):
      _observation_order(order),
      _seq(move(seq)){}
   };
   int _next_obs_order;
   bool _keep_seq;
   typedef unordered_map<string, SequenceInfo> id2SeqInfo;
   id2SeqInfo _by_id;

public:

   RefSeqTable(bool keep_seq);
   bool insertIntoTable(const string name, unique_ptr<DNABitSet> seq);
   unique_ptr<DNABitSet> get_seq(const string name);
   int observation_order(const string name) const;
};


class HitFactory
{
private:
   static const size_t kHitBufMaxSize = 10 * 1024;
   ReadTable& _reads_table;
   RefSeqTable& _ref_table;
   char _hit_buf[kHitBufMaxSize];

public:
   HitFactory(ReadTable &reads_table, RefSeqTable &ref_table);
   HitFactory(HitFactory &rhs) = delete; //non-copible class.
   HitFactory& operator=(HitFactory &rhs) = delete; // non-copible class
   HitFactory(HitFactory &&rhs) = default;
   HitFactory& operator=(HitFactory &&rhs) = default;
   virtual ~HitFactory() = default;
   virtual bool recordsRemain() const = 0;
   virtual bool nextRecord(const char* &buf, size_t& buf_size) = 0;
   virtual bool getHitFromBuf(const char* bwt_buf, ReadHit& bh)=0;
   virtual RefSeqTable& ref_table() { return _ref_table; }
   virtual ReadTable& reads_table(){return _reads_table;}
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
   bool recordsRemain() const;
   void markCurrPos();
   void reset();
   void undo_hit();
   bool nextRecord(const char* &buf, size_t& buf_size);
   bool getHitFromBuf(const char* bwt_buf, ReadHit& bh);
};

typedef unique_ptr<const ReadHit> ReadHitPtr;

//PairedHit solely own the ReadHit objects by
//using unique_ptr
class PairedHit{
   ReadHitPtr _left_read ;
   ReadHitPtr _right_read;
   double _collapse_mass = 0.0;

public:
   PairedHit() = default;
   PairedHit(ReadHitPtr leftRead, ReadHitPtr rightRead);
   ReadHitPtr left_read();
   void set_left_read(ReadHitPtr lr);
   ReadHitPtr right_read();
   void set_right_read(ReadHitPtr rr);
   bool is_paired() const;
   string ref_seq_name() const;
   int left_pos() const;
   int right_pos() const;
   uint edit_dist() const;
   bool contains_splice() const;
   ReadID read_id() const;

   bool paried_hit_lt(const PairedHit &rhs) const;
};


#endif /* STRAWB_READ_H_ */
