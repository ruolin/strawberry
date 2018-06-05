/*
 * This file is the interface that reads BAM/SAM files.
 * Utilize some functions from Cufflinks:
 * https://github.com/cole-trapnell-lab/cufflinks/blob/master/src/hits.h
 */

#ifndef READ_HPP
#define READ_HPP
#include<vector>
#include<memory>
#include<string>
#include<unordered_map>
#include<cassert>
#include <iostream>
#include "common.h"
#include "sam/sam.h"
//#include "kmer.h"

//using namespace std;
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

inline std::ostream& operator<<(std::ostream& os, const CigarOp& co){
   switch(co._type){
      case BAM_CMATCH:
         os<<co._length<<"M";
         break;
      case INS:
         os<<co._length<<"I";
         break;
      case DEL:
         os<<co._length<<"D";
         break;
      case REF_SKIP:
         os<<co._length<<"N";
         break;
      case SOFT_CLIP:
         os<<co._length<<"S";
         break;
      case HARD_CLIP:
         os<<co._length<<"H";
         break;
      default:
         break;
   }
   return os;
}

class ReadHit{
private:

   ReadID _read_id;
   std::string _read_name;
   GenomicInterval _iv;
   std::vector<CigarOp> _cigar;
   RefID _partner_ref_id;
   uint _partner_pos;
   int _num_mismatch = -1;
   int _num_hits = 1;
   int _trans_left=0; // position relate to transcriptom
   uint32_t _sam_flag = 0; //bitwise FLAG
   double _read_mass = 0.0;
public:
   std::string _seq;
   ReadHit() = default;
   ReadHit( ReadID readID,
            std::string read_name,
         GenomicInterval iv,
         const std::vector<CigarOp> & cigar,
         RefID partnerRef,
         int partnerPos,
         int numMismatch,
         int numHit,
         uint32_t samFlag,
         double mass,
         char* seq
         );
//   ~ReadHit(){
//      if(_seq != NULL){
//         free(_seq);
//         _seq = NULL;
//      }
//   }
   uint read_len() const;
   uint intron_lens() const; /*Return intron len, otherwise return 0*/
   std::vector<std::pair<uint,uint>> intron_coords() const;
   ReadID read_id() const;
   bool contains_splice() const;
   GenomicInterval interval() const;
   RefID partner_ref_id() const;
   const std::vector<CigarOp>& cigar() const;
   uint partner_pos() const;
   RefID ref_id() const; // chromosome id or scaffold id containing the read
   int num_mismatch() const;
   int numHits() const;
   bool is_singleton() const;
   std::string read_name() const {
      return _read_name;
   }
   uint left() const;
   uint right() const;
   Strand_t strand() const;
   double mass() const;
   double raw_mass() const;
   void mass(double m);
   bool is_first() const;
   bool is_second() const;

   bool reverseCompl() const;
   bool operator==(const ReadHit& rhs) const; // not considering read orientation
   bool operator!=(const ReadHit& rhs) const;
   bool operator<(const ReadHit& rhs) const;
   //char* read_seq() const;
   //std::vector<CigarOp> cigars() const;
};

// Now it is used only to convert read name to read id
//
class ReadTable
{
 public:
	 
   // This function should NEVER return zero
   ReadID get_id(const std::string& name);
   uint _read_len_abs = 0;
   std::vector<int> _frag_dist;
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

   int _total_reads;
   std::vector<double> _emp_dist;
   double _mean;
   double _sd;
   bool _use_emp;
   int _start_offset;
   int _end_offset;
   InsertSize();
   InsertSize(double mean, double sd);
   InsertSize(const std::vector<int> frag_lens);
   double emp_dist_pdf(uint insert_size) const;
   bool empty() const;
   //double truncated_normal_pdf(uint insert_size) const;
};


class RefSeqTable
{

private:
//_id is the observation order.
//_id start from 0 which is used as the index in std::vector<GffSeqData>;
   std::string _seq;
   bool _keep_seq;
   std::vector<std::string> _id_2_real_name;
   std::vector<std::string> _id2name;
   std::unordered_map<std::string, int> _name2id;
public:
   RefSeqTable(bool keep_seq) : _keep_seq(keep_seq){}
   int get_id(std::string& name);
   int set_id(std::string& name);
   int size() const{
      return _name2id.size();
   }
   const std::string ref_real_name(int id) const;
};


class HitFactory
{
private:
   virtual platform_t str2platform(const std::string pl_str);
protected:
   static const unsigned MAX_HEADER_LEN = 4 * 1024 * 1024; // 4 MB
   static const size_t kHitBufMaxSize = 10 * 1024;
   char _hit_buf[kHitBufMaxSize];
   int _num_seq_header_recs = 0;
   std::string _hit_file_name;
   AssayProperties _assay_props;
public:
   ReadTable& _reads_table;
   RefSeqTable& _ref_table;
   HitFactory(ReadTable &reads_table, RefSeqTable &ref_table, std::string hit_file_name);
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
   virtual bool parse_header_line(const std::string& hline);
   virtual bool inspect_header() = 0;
   virtual void reset() = 0;
   virtual std::string sample_name() const {
      return _hit_file_name;
   }
   virtual void return2Pos(int64_t pos) = 0;
   virtual int64_t getCurrPos() = 0;
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
   BAMHitFactory(const std::string& bam_file_name,
                 ReadTable& read_table,
                 RefSeqTable &ref_table);
   ~BAMHitFactory();
   bool recordsRemain() const;
   void markCurrPos();
   void reset();
   void undo_hit();
   bool nextRecord(const char* &buf, size_t& buf_size);
   bool getHitFromBuf(const char* bwt_buf, ReadHit& bh);
   bool inspect_header();
   int64_t getCurrPos();
   void return2Pos(int64_t pos);
};

typedef std::shared_ptr<ReadHit> ReadHitPtr;


class PairedHit{
   double _collapse_mass = 0.0;
   double _mass = 0.0;
public:
   ReadHitPtr _right_read ;
   ReadHitPtr _left_read ;
   //std::vector<Kmer> _left_kmers;
   //std::vector<Kmer> _right_kmers;

   PairedHit() = default;
   PairedHit(ReadHitPtr leftRead, ReadHitPtr rightRead);
   PairedHit& operator=(const PairedHit& rhs){
      _left_read = rhs._left_read;
      _right_read = rhs._right_read;
      _collapse_mass = rhs._collapse_mass;
      return *this;
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
   double raw_mass() const;
   bool operator==(const PairedHit& rhs) const;
   bool operator!=(const PairedHit& rhs) const;
   bool operator<(const PairedHit& rhs) const;
   bool paried_hit_lt(const PairedHit &rhs) const;
   void add_2_collapse_mass(double add);
   double collapse_mass() const;
   void init_raw_mass();
   void weighted_mass(double m);
   double weighted_mass() const;
   //void set_kmers(int num_kmers);
};

void mean_and_sd_insert_size(const std::vector<int> & vec, double & mean, double &sd);
#endif /* READ_HPP */
