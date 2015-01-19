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
#include<GBase.h>
#include<unordered_map>

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
   char _strand; // strand plus/minus known for spliced read, otherwise unknown.
   int _num_mismatch = -1;
   int _num_hit = 1;
   int _trans_left=0; // position relate to transcriptom
   uint32_t _sam_flag = 0;
   float _read_mass = 0.0;

public:
   ReadHit();
   ReadHit( ReadID readID,
         GenomicInterval iv,
         const vector<CigarOp> & cigar,
         char strand,
         string partnerRef,
         int partnerPos,
         int numMismatch,
         int numHit,
         uint32_t samFlag,
         float mass,
         );

   int read_len() const;
   bool containsSplice() const;
   string partner_ref() const;
   int partner_pos() const;
   string ref() const; // chromosome or scaffold containing the read
   bool is_singleton() const;
   int left() const;
   int right() const;
   char strand() const;
};

typedef shared_ptr<const ReadHit> ReadHitPtr;

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

class RefSequenceTable
{
private:

   struct SequenceInfo{
      uint32_t _observation_order;
      unique_ptr<string> _seq;
      SequenceInfo(uint32_t order,
                   string name,
                   unique_ptr<DNABitSet> seq):
      _observation_order(order),
      _seq(seq){}
   };

   bool _keep_seq;
   typedef unordered_map<string, SequenceInfo> id2SeqInfo;
   id2SeqInfo _by_id;

public:

   RefSequenceTable(bool keep_seq);
   insertIntoTable(const string  name, unique_ptr<string> seq);
   const string get_name(RefID ID) const;
   unique_ptr<string> get_seq(RefID ID) const;
   //const unique_ptr<SequenceInfo> get_info(RefID ID) const;
   int observation_order(RefID ID) const;
};


class HitFactory{
private:
   static const size_t kHitBufMaxSize = 10 * 1024;
   ReadTable &_reads_table;
   RefSequenceTable & _ref_table;
   char _hit_buf[kHitBufMaxSize];

public:
   HitFactory(ReadTable &reads_table, RefSequenceTable ref_table);
   HitFactory(HitFactory &rhs) = delete; //non-copible class.
   HitFactory& operator=(HitFactory &rhs) = delete; // non-copible class
   HitFactory(HitFactory &&rhs);
   HitFactory& operator=(HitFactory &&rhs);
   virtual ~HitFactory() = default;
   virtual bool recordsRemain() const = 0;
   virtual void markCurrPos();
   virtual bool nextRecord(const char* &buf, size_t& buf_size) = 0;
   virtual bool getHitFromBuf(const char* bwt_buf, ReadHit& bh)=0;
};

class BAMHitFactory {
   int64_t _curr_pos;
   samfile_t* _hit_file;
   bam1_t _next_hit;
   bool _eof_encountered;
   int64_t _beginning;
public:
   BAMHitFactory(const string& bam_file_name);
   ~BAMHitFactory();
   bool recordsRemain() const;
   void markCurrPos();
   bool nextRecord(const char* &buf, size_t& buf_size);
   bool getHitFromBuf(const char* bwt_buf, ReadHit& bh);
   void reset();
};

class PairedHit{
   ReadHitPtr _left_read ;
   ReadHitPtr _right_read;
public:
   PairedHit(ReadHitPtr leftRead, ReadHitPtr rightRead);
   ReadHitPtr getLeftRead() const;
   void setLeftRead(ReadHitPtr lr);
   ReadHitPtr getRightRead() const;
   void setRightRead(ReadHitPtr lr);
   int getLeftPos() const;
   int getRightPos() const;
   int genomicInnerDist() const;
   uint edit_dist() const;
   bool paried_hit_lt(const PairedHit &rhs) const;
};


#endif /* STRAWB_READ_H_ */
