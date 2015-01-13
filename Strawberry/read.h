/*
 * read.h
 *
 *  Created on: Nov 3, 2014
 *      Author: RUOLIN
 */

#ifndef STRAWB_READ_H_
#define STRAWB_READ_H_
#include<stdio.h>
#include<string>
#include<bam/sam.h>
#include <cassert>
#include<vector>
#include<memory>


using namespace std;
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
   CigarOp(CigarOpCode o, uint32_t l) : opcode(o), length(l) {}
   CigarOpCode opcode : 3;
   uint32_t length : 29;

   bool operator==(const CigarOp& rhs) const { return opcode == rhs.opcode && length == rhs.length; }

};

typedef uint64_t ReadID;
class ReadHit{
public:
   bool _left_most; // left most reads in transcriptional direction
   bool _paired;
   int _start; // left of read span
   vector<CigarOp> _cigar;
   int _num_mismatch;
   int _num_hit;

   string _ref_id;
   ReadID _read_id;
   int _trans_left=0; // position relate to transcriptom
   char _strand;
   string _partner_ref;
   int _partner_pos;
   uint32_t _sam_flag;
   ReadHit()=default;
   ReadHit( ReadID readID,
         string refID,
         int pos,
         const vector<CigarOp> & cigar,
         char strand,
         string partnerRef,
         int partnerPos,
         int numMismatch,
         int numHit,
         uint32_t samFlag,
         bool paired
         );
   ReadHit(const ReadHit& other);
   ~ReadHit(){};
   int read_len() const;
   int getEnd() const;
   void setStrandness(char gene_strand);
   int getStart() const;
};

typedef shared_ptr<const ReadHit> ReadHitPtr;

inline ReadID string2Hash(const string &name){
   const char* __s = name.c_str();
   uint64_t hash = 0xcbf29ce484222325ull;
   for( ; *__s; ++__s){
      hash *= 1099511628211ull;
      hash ^= *__s;
   }
   assert(hash);
   return hash;
}

/*
class HitFactory{
   friend ReadHit;
private:
   static const size_t _hit_buf_max_sz = 10 * 1024;
   char _hit_buf[_hit_buf_max_sz];
   vector<ReadHit> _hits={};
public:
   HitFactory()= default;
   ~HitFactory()
   virtual bool recordsRemain() const = 0;
   virtual void markCurrPos();
   virtual bool nextRecord(const char* &buf, size_t& buf_size) = 0;
   virtual bool getHitFromBuf(const char* bwt_buf, ReadHit& bh)=0;
};
*/
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
