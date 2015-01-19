/*
 * read.cpp
 *
 *  Created on: Nov 3, 2014
 *      Author: RUOLIN LIU
 */

#include "read.h"
#include <algorithm>
#include<assert.h>
#include <stdexcept>

ReadHit::ReadHit():
   _read_id(0),
   _partner_pos(0),
   _strand('.') {}


ReadHit::ReadHit(
   ReadID readID,
   GenomicInterval iv,
   const vector<CigarOp> & cigar,
   char strand,
   string partnerRef,
   int partnerPos,
   int numMismatch,
   int numHit,
   uint32_t samFlag,
   float mass ):
      _read_id(readID),
      _iv(iv),
      _cigar(cigar),
      _strand(strand),
      _partner_ref(partnerRef),
      _partner_pos(partnerPos),
      _num_mismatch(numMismatch),
      _num_hit(numHit),
      _sam_flag(samFlag)
{}



int ReadHit::read_len() const {return _iv.len();}

bool ReadHit::containsSplice()const{
      for (size_t i = 0; i < _cigar.size(); ++i){

            if (_cigar[i].opcode == REF_SKIP)
               return true;
      }
      return false;
}

string ReadHit::ref() const {return _iv.chrom();}

string ReadHit::partner_ref() const { return _partner_ref;}

int ReadHit::partner_pos() const { return _partner_pos;}

int ReadHit::right() const {return _iv.right();}

char ReadHit::strand() const {return _iv.strand();}

int ReadHit::left() const { return _iv.left();}

bool ReadHit::is_singleton() const
{
   return (partner_ref() == 0 ||
         partner_ref() != ref() ||
         abs(partner_pos() - left()) > max_partner_dist);
}




ReadID ReadTable::get_id(const string& name)
{
   uint64_t id = hashString(name.c_str());
   assert(id);
   return id;
}

RefSequenceTable::RefSequenceTable(bool keep_seq):
   _keep_seq(keep_seq),
   _next_id(1){}

RefID RefSequenceTable::init_id(const string& name, unique_ptr<string> seq){
   uint64_t id = hashString(name.c_str());
   pair<id2SeqInfo::iterator, bool> res = _by_id.insert(make_pair(id,
         SequenceInfo(_next_id, name, nullptr) )  );
   if(res.second){
      if(_keep_seq) res.first->second._seq = seq;
      ++_next_id;
   }
   assert (id);
   return id;
}

const string RefSequenceTable::get_name(RefID ID) const{
   auto it = _by_id.find(ID);
   if( it != _by_id.end()) return it->second._name;
   else GError("ID %d is not in the Reference Sequence Table\n", ID);
}

unique_ptr<string> RefSequenceTable::get_seq(RefID ID) const{
   auto it = _by_id.find(ID);
   if( it != _by_id.end()) return it->second._seq;
   else GError("ID %d is not in the Reference Sequence Table\n", ID);
}

int RefSequenceTable::observation_order(RefID ID) const{
   auto it = _by_id.find(ID);
   if( it != _by_id.end()) return it->second._observation_order;
   else GError("ID %d is not in the Reference Sequence Table\n", ID);
}
//const unique_ptr<RefSequenceTable::SequenceInfo> RefSequenceTable::get_info(RefID ID) const{
//   auto it = _by_id.find(ID);
//   if( it != _by_id.end()) return &(it->second);
//   else GError("ID %d is not in the Reference Sequence Table\n", ID);
//}

HitFactory::HitFactory(ReadTable &reads_table, RefSequenceTable ref_table):
   _reads_table(reads_table),_ref_table(ref_table){}

HitFactory::HitFactory(HitFactory &&rhs):
   _reads_table( move(rhs._reads_table)), _ref_table(move(rhs._ref_table)){}

HitFactory& HitFactory::operator =(HitFactory &&rhs){
   if(this != &rhs){
      _reads_table = move(rhs._reads_table);
      _ref_table = move(rhs._ref_table);
   }
   return *this;
}
BAMHitFactory::BAMHitFactory(const string& bam_file_name){
   _hit_file = samopen(bam_file_name.c_str(), "rb", 0);
   memset(&_next_hit, 0, sizeof(_next_hit));
   if(_hit_file == NULL || _hit_file->header == NULL){
      throw std::runtime_error("Fail to open BAM file");
   }
   _beginning = bgzf_tell(_hit_file->x.bam);
   _eof_encountered = false;

}

BAMHitFactory::~BAMHitFactory(){
   if(_hit_file) samclose(_hit_file);
}

void BAMHitFactory::reset()
   {
      if (_hit_file && _hit_file->x.bam)
      {
         bgzf_seek(_hit_file->x.bam, _beginning, SEEK_SET);
            _eof_encountered = false;
      }
   }
void BAMHitFactory::markCurrPos(){
   _curr_pos = bgzf_tell(_hit_file->x.bam);
}

bool BAMHitFactory::nextRecord(const char* &buf, size_t& buf_size){
   if(_next_hit.data){
      free(_next_hit.data);
      _next_hit.data=NULL;
   }
   if (recordsRemain() == false)
      return false;
   markCurrPos();
   memset(&_next_hit,0,sizeof(_next_hit));
   int bytes_read = samread(_hit_file, &_next_hit);
   if (bytes_read < 0){
      _eof_encountered = true;
      return false;
   }
   buf = (const char*)& _next_hit;
   buf_size = bytes_read;
}

bool BAMHitFactory::recordsRemain() const{
   return !_eof_encountered;
}

bool BAMHitFactory::getHitFromBuf(const char* orig_bwt_buf, ReadHit &bh){
   const bam1_t* hit_buf = (const bam1_t*)orig_bwt_buf;
   uint32_t sam_flag = hit_buf->core.flag;
   int pos = hit_buf->core.pos;
   int mate_pos = hit_buf->core.mpos;
   int target_id = hit_buf->core.tid;
   int mate_target_id = hit_buf->core.mtid;
   ReadID readid = string2Hash(bam1_qname(hit_buf));
   vector<CigarOp> cigar ={};
   int num_hits = 1;
   string text_name = _hit_file->header->target_name[target_id];
   if( (sam_flag & 0x4) || target_id) return true;  // unmapped
   bool paired = sam_flag & BAM_FPAIRED;
   bool antisense_aln = sam_flag & 0x10;
   if(target_id >= _hit_file->header->n_targets){
      fprintf(stderr, "BAM error: file contains hits to sequences not in BAM file header");
      return false;
   }
   for (int i=0; i<hit_buf->core.n_cigar; ++i){
      int length = bam1_cigar(hit_buf)[i] >> BAM_CIGAR_SHIFT;
      if(length <= 0){
         fprintf(stderr, "BAM error: CIGAR op has zero length (%s)\n",bam1_qname(hit_buf));
         return false;
      }
      CigarOpCode opcode;
      switch(bam1_cigar(hit_buf)[i] & BAM_CIGAR_MASK){
      case BAM_CMATCH: opcode = MATCH; break;
      case BAM_CINS: opcode = INS; break;
      case BAM_CDEL: opcode = DEL; break;
      case BAM_CSOFT_CLIP: opcode = SOFT_CLIP; break;
      case BAM_CHARD_CLIP: opcode = HARD_CLIP; break;
      case BAM_CPAD: opcode = HARD_CLIP; break;
      case BAM_CREF_SKIP: opcode = REF_SKIP;
            // if(length > (int)max_intron_length) return false; // for future
            break;
      default: return false;
      }
      if (opcode != HARD_CLIP)
         cigar.push_back(CigarOp(opcode, length));
   }

   string mate_ref;
   if(mate_target_id >=0){
      if(mate_target_id == target_id){
         mate_ref = _hit_file->header->target_name[mate_target_id];
      }
      else return false;
   }
   else mate_pos = 0;

   char source_strand = '.';
   unsigned char num_mismatches = 0;
   if(antisense_aln){
      source_strand = '-';
   }
   else
      source_strand = '+';

   uint8_t* ptr = bam_aux_get(hit_buf, "NM");
   if(ptr){
      num_mismatches = bam_aux2i(ptr);
   }

   ptr = bam_aux_get(hit_buf, "NH");
   if(ptr){
      num_hits = bam_aux2i(ptr);
   }

   /*
    * havn't use mass property
   double mass = 1.0;
       ptr = bam_aux_get(hit_buf, "ZF");
      if (ptr)
      {
         mass = bam_aux2i(ptr);
           if (mass <= 0.0)
               mass = 1.0;
      }
   */

   bh = ReadHit(readid,
         text_name,
            pos,
            cigar,
            source_strand,
            mate_ref,
            mate_pos,
            num_mismatches,
            num_hits,
            sam_flag,
            paired
            );
}


PairedHit::PairedHit(ReadHitPtr leftRead, ReadHitPtr rightRead):
            _left_read(leftRead), _right_read(rightRead){}

ReadHitPtr PairedHit::getLeftRead() const {return _left_read;}

void PairedHit::setLeftRead(ReadHitPtr lr) {_left_read = lr;}

ReadHitPtr PairedHit::getRightRead() const {return _right_read;}

int PairedHit::getLeftPos() const{
   if(_right_read && _left_read){
      return min(_right_read->getStart(), _left_read->getStart());
   }
   else if (_left_read)
      return _left_read->getStart();
   else if (_right_read)
      return _right_read->getStart();
   else
      return -1;
}

int PairedHit::getRightPos() const{
   if(_right_read && _left_read){
      return max(_right_read->getEnd(), _left_read->getEnd());
   }
   else if(_right_read)
      return _right_read->getEnd();
   else if(_left_read)
      return _left_read->getEnd();
   else
      return -1;
}

int PairedHit::genomicInnerDist() const{
   if(_left_read && _right_read)
      return _right_read->getStart() - _left_read->getEnd();
   return -1;
}

uint PairedHit::edit_dist() const{
   uint edits = 0;
   if(_left_read)
      edits += _left_read->_num_mismatch;
   if(_right_read)
      edits += _right_read->_num_mismatch;
   return edits;
}

bool PairedHit::paried_hit_lt (const PairedHit &rhs) const{
   if( getLeftPos() != rhs.getLeftPos() )
      return getLeftPos() < rhs.getLeftPos();
   if( getRightPos() != rhs.getRightPos() )
      return getRightPos() > rhs.getRightPos();
   // proceed only when same start and end;
   if ((_left_read == NULL) != (rhs.getLeftRead()==NULL) )
      return (_left_read == NULL) < (rhs.getLeftRead() == NULL);
   if ((_right_read ==NULL) != (rhs.getRightRead()==NULL) )
      return (_right_read == NULL) < (rhs.getRightRead() == NULL);
   assert ((_left_read == NULL) == (rhs.getLeftRead()==NULL) );
   assert ((_right_read ==NULL) == (rhs.getRightRead()==NULL) );
}


