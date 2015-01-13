/*
 * read.cpp
 *
 *  Created on: Nov 3, 2014
 *      Author: RUOLIN LIU
 */

#include "read.h"
#include <algorithm>
#include <stdexcept>
ReadHit::ReadHit(
   ReadID readID,
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
   ) :_ref_id(refID),
      _read_id(readID),
      _start(pos),
      _cigar(cigar),
      _strand(strand),
      _partner_ref(partnerRef),
      _partner_pos(partnerPos),
      _num_mismatch(numMismatch),
      _num_hit(numHit),
      _paired(paired),
      _sam_flag(samFlag)
{}


ReadHit::ReadHit(const ReadHit& other):
   _ref_id(other._ref_id),
   _read_id(other._read_id),
   _start(other._start),
   _cigar(other._cigar),
   _trans_left(other._trans_left),
   _strand(other._strand),
   _partner_ref(other._partner_ref),
   _partner_pos(other._partner_pos),
   _num_mismatch(other._num_mismatch),
   _num_hit(other._num_hit),
   _paired(other._paired),
   _sam_flag(other._sam_flag)
   {}

int ReadHit::read_len() const{
      int len = 0;
      for (size_t i = 0; i < _cigar.size(); ++i)
      {
         const CigarOp& op = _cigar[i];
         switch(op.opcode)
         {
            case MATCH:
            case INS:
            case SOFT_CLIP:
               len += op.length;
               break;
            default:
               break;
         }
      }

      return len;
   }

int ReadHit::getStart() const{
   return _start;
}

int ReadHit::getEnd() const{
      int r = _start;
      for (size_t i = 0; i < _cigar.size(); ++i)
      {
         const CigarOp& op = _cigar[i];

         switch(op.opcode)
         {
            case MATCH:
            case REF_SKIP:
                case SOFT_CLIP:
            case DEL:
               r += op.length;
               break;
                case INS:
                case HARD_CLIP:
            default:
               break;
         }
      }
      return r;
}

void ReadHit::setStrandness(char gene_strand){
   if (gene_strand == '+'){
      if(_strand == '+') _left_most = true;
      else if (_strand ==  '-' ) _left_most = false;
   }
   else if(gene_strand == '-'){
      if(_strand == '+') _left_most = false;
      else if(_strand == '-') _left_most = true;
   }
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
