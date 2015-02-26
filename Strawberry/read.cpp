/*
 * read.cpp
 *
 *  Created on: Nov 3, 2014
 *      Author: RUOLIN LIU
 */

#include <algorithm>
#include<assert.h>
#include <stdexcept>
#include <string.h>
#include "read.h"

static const int kMaxIntronLength = 60000;
ReadHit::ReadHit(
   ReadID readID,
   GenomicInterval iv,
   const vector<CigarOp> & cigar,
   RefID partnerRef,
   int partnerPos,
   int numMismatch,
   int numHit,
   uint32_t samFlag,
   float mass ):
      _read_id(readID),
      _iv(iv),
      _cigar(cigar),
      _partner_ref_id(partnerRef),
      _partner_pos(partnerPos),
      _num_mismatch(numMismatch),
      _num_hits(numHit),
      _sam_flag(samFlag)
{}



uint ReadHit::read_len() const { return _iv.len();}

double ReadHit::mass() const{
   if(is_singleton()){
      return 1.0/_num_hits;
   }else{
      return 0.5/_num_hits;
   }
}

bool ReadHit::contains_splice()const{
      for (size_t i = 0; i < _cigar.size(); ++i){

            if (_cigar[i].opcode == REF_SKIP)
               return true;
      }
      return false;
}

bool ReadHit::is_first() const
{
   return _sam_flag & 0x40;
}

ReadID ReadHit::read_id() const {return _read_id;}

RefID ReadHit::ref_id() const {return _iv.seq_id();}

RefID ReadHit::partner_ref_id() const { return _partner_ref_id;}

uint ReadHit::partner_pos() const { return _partner_pos;}

uint ReadHit::right() const {return _iv.right();}

GenomicInterval ReadHit::interval() const { return _iv;}

char ReadHit::strand() const {return _iv.strand();}

//vector<CigarOp> ReadHit::cigars() const {
//   return _cigar;
//}

uint ReadHit::left() const { return _iv.left();}

int ReadHit::numHits() const { return _num_hits;}
int ReadHit::num_mismatch() const { return _num_mismatch;}

bool ReadHit::is_singleton() const
{
   return (partner_pos() == 0 ||
         partner_ref_id() == -1 ||
         partner_ref_id() != ref_id());
}

bool ReadHit::reverseCompl() const
{
   return _sam_flag  & 0x10;
}

// not considering read orientation
bool ReadHit::operator==(const ReadHit& rhs) const
{
   assert(!_cigar.empty() && !rhs._cigar.empty());
   return (_iv == rhs._iv && rhs._cigar == _cigar);
}



ReadID ReadTable::get_id(const string& name)
{
   uint64_t id = hashString(name.c_str());
   assert(id);
   return id;
}


//const unique_ptr<RefSeqTable::SequenceInfo> RefSeqTable::get_info(RefID ID) const{
//   auto it = _by_id.find(ID);
//   if( it != _by_id.end()) return &(it->second);
//   else GError("ID %d is not in the Reference Sequence Table\n", ID);
//}

HitFactory::HitFactory(ReadTable &reads_table, RefSeqTable &ref_table):
   _reads_table(reads_table),_ref_table(ref_table){}

platform_t HitFactory::str2platform(const string str)
{
    if (str == "SOLiD")
    {
        return SOLID;
    }
    else if (str == "Illumina")
    {
        return ILLUMINA;
    }
    else
    {
        return UNKNOWN_PLATFORM;
    }
}

bool HitFactory::parse_header_line(const string& hline){
//#ifdef DEBUG
//   cout<<hline<<endl;
//#endif
   vector<string> cols;
   split(hline, "\t", cols);
   if(cols[0] == "@SQ"){
      ++_num_seq_header_recs;
      for(auto &i: cols){
         vector<string> fields;
         split(i, ":", fields);
         if(fields[0] == "SN"){
            str2lower(fields[1]);
            RefID _ID = _ref_table.get_id(fields[1]);
            if(_ID != _num_seq_header_recs-1)
               SMessage("Error: sort order of reads in BAM not consistent.\n");
               return false;
         }
      }
   }

   if(cols[0] == '@RG'){
      for(auto &i: cols){
         vector<string> fields;
         split(i, ":", fields);
         if(fields[0] == "PL"){
            platform_t p = str2platform(fields[1]);
            assert(_read_group_props._platform == UNKNOWN_PLATFORM);
            _read_group_props._platform = p;
         }
      }
   }
   return true;
}

const ReadGroupProperties&  HitFactory::read_group_properties()
{
   return _read_group_props;
}

BAMHitFactory::BAMHitFactory(const string& bam_file_name,
                            ReadTable &read_table,
                            RefSeqTable &ref_table) throw(runtime_error):
                 HitFactory(read_table, ref_table)
{
   _hit_file = samopen(bam_file_name.c_str(), "rb", 0);
   memset(&_next_hit, 0, sizeof(_next_hit));
   if(_hit_file == NULL || _hit_file->header == NULL){
      throw runtime_error("Fail to open BAM file");
   }

   _beginning = bgzf_tell(_hit_file->x.bam);
   _curr_pos = _beginning;
   _eof_encountered = false;

}

BAMHitFactory::~BAMHitFactory()
{
   if (_hit_file){
      samclose(_hit_file);
      free(_next_hit.data);
   }
}

bool BAMHitFactory::inspect_header()
{
   bam_header_t* header = _hit_file->header;
   if(header == NULL || header->l_text == 0){
      SMessage("Warning: No BAM header\n");
      return false;
   }

   if(header->l_text >= MAX_HEADER_LEN ){
      SMessage("Warning: BAM header too large\n");
      return false;
   }

   if(header->text != NULL){
      char* h_text = strdup(header->text);
      char* pBuf = h_text;
      while( pBuf - h_text < header->l_text){
         char *nl = strchr(pBuf,'\n');
         if (nl){
            *nl = 0;
            parse_header_line(pBuf);
            pBuf = ++nl;
         }
         else{
            pBuf = h_text + header->l_text;
         }
      }

      free(h_text);
   }
   return true;
}

void BAMHitFactory::reset()
{
   if (_hit_file && _hit_file->x.bam)
   {
      bgzf_seek(_hit_file->x.bam, _beginning, SEEK_SET);
      _eof_encountered = false;
   }
}

void BAMHitFactory::markCurrPos()
{
   _curr_pos = bgzf_tell(_hit_file->x.bam);
}

bool BAMHitFactory::recordsRemain() const{
   return !_eof_encountered;
}

void BAMHitFactory::undo_hit(){
   bgzf_seek(_hit_file->x.bam, _curr_pos, SEEK_SET);
}

bool BAMHitFactory::nextRecord(const char* &buf, size_t& buf_size)
{
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

   return true;
}

bool BAMHitFactory::getHitFromBuf(const char* orig_bwt_buf, ReadHit &bh){
   const bam1_t* hit_buf = (const bam1_t*)orig_bwt_buf;
   uint32_t sam_flag = hit_buf->core.flag;
   uint pos = hit_buf->core.pos + 1; // BAM file index starts at 0
   uint mate_pos = hit_buf->core.mpos + 1; // BAM file index starts at 0
   int target_id = hit_buf->core.tid;
   int mate_target_id = hit_buf->core.mtid;
   string mate_text_name;

   if(mate_target_id < 0){
      mate_text_name = "*";
   }
   else{
      mate_text_name = _hit_file->header->target_name[mate_target_id];
      str2lower(mate_text_name);
   }
   RefID parterner_ref_id = ref_table().get_id(mate_text_name);

   int read_len = 0;
   ReadID readid = HitFactory::reads_table().get_id(bam1_qname(hit_buf));
   vector<CigarOp> cigar;
   bool is_spliced_alignment = false;
   int num_hits = 1;
   if( (sam_flag & 0x4) || target_id < 0 ){
      bh = ReadHit(readid,
                   GenomicInterval(),
                   cigar,
                   parterner_ref_id,
                   0,
                   0,
                   num_hits,
                   sam_flag,
                   0.0
                   );
      return true;
   }


   string text_name = _hit_file->header->target_name[target_id];
   str2lower(text_name);
   RefID ref_id = ref_table().get_id(text_name);

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
      case BAM_CMATCH: opcode = MATCH;
         read_len += length;
         break;
      case BAM_CINS: opcode = INS;
         read_len += length;
         break;
      case BAM_CDEL: opcode = DEL; break;
      case BAM_CSOFT_CLIP: opcode = SOFT_CLIP;
         read_len += length;
         break;
      case BAM_CHARD_CLIP: opcode = HARD_CLIP; break;
      case BAM_CPAD: opcode = HARD_CLIP; break;
      case BAM_CREF_SKIP:
         opcode = REF_SKIP;
         is_spliced_alignment = true;
         if(length > (int) kMaxIntronLength){
            SMessage("length %d is larger than max intron length %d \n", length, kMaxIntronLength);
            return false;
         }
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
   } else{
      mate_pos = 0;
   }

   char source_strand = GenomicInterval::kStrandUnknown;
   unsigned char num_mismatches = 0;

   uint8_t* ptr = bam_aux_get(hit_buf, "XS");
   if(ptr){
      char src_strand_char = bam_aux2A(ptr);
      source_strand = src_strand_char;
   }

   ptr = bam_aux_get(hit_buf, "NM");
   if(ptr){
      num_mismatches = bam_aux2i(ptr);
   }

   ptr = bam_aux_get(hit_buf, "NH");
   if(ptr){
      num_hits = bam_aux2i(ptr);
   }

   double mass = 1.0;
   ptr = bam_aux_get(hit_buf, "ZF");
   if (ptr)
   {
      mass = bam_aux2i(ptr);
        if (mass <= 0.0)
            mass = 1.0;
   }

   if(is_spliced_alignment){
      if(source_strand == GenomicInterval::kStrandUnknown)
         fprintf(stderr, "BAM record error: found spliced alignment without XS attribute\n");
   }

   bh = ReadHit(
               readid,
               GenomicInterval(ref_id, pos, pos+read_len-1, source_strand),
               cigar,
               parterner_ref_id,
               mate_pos,
               num_mismatches,
               num_hits,
               sam_flag,
               mass
               );
   return true;
}


PairedHit::PairedHit(ReadHitPtr leftRead, ReadHitPtr rightRead):
            _left_read(move(leftRead)), _right_read(move(rightRead)){}

const ReadHitPtr PairedHit::left_read() const {return move(_left_read);}

char PairedHit::strand() const {
   if(_right_read && _left_read){
      assert(_right_read->strand() == _left_read->strand() ||
            _right_read->strand() == GenomicInterval::kStrandUnknown ||
            _left_read->strand() == GenomicInterval::kStrandUnknown);
         return _left_read->strand();
   } else if (_left_read)
      return _left_read->strand();
   else if (_right_read)
      return _right_read->strand();
   else
      return 0;
}

void PairedHit::set_left_read(ReadHitPtr lr)
{
   _left_read = move(lr);
}

const ReadHitPtr PairedHit::right_read() const {return move(_right_read);}

void PairedHit::set_right_read(ReadHitPtr rr)
{
   _right_read = move(rr);
}

uint PairedHit::left_pos() const{
   if(_right_read && _left_read){
      return min(_right_read->left(), _left_read->left());
   }
   else if (_left_read)
      return _left_read->left();
   else if (_right_read)
      return _right_read->left();
   else
      return -1;
}

uint PairedHit::right_pos() const{
   if(_right_read && _left_read){
      return max(_right_read->right(), _left_read->right());
   }
   else if(_right_read)
      return _right_read->right();
   else if(_left_read)
      return _left_read->right();
   else
      return -1;
}

bool PairedHit::is_paired() const
{
   return _left_read && _right_read;
}



uint PairedHit::edit_dist() const{
   uint edits = 0;
   if(_left_read)
      edits += _left_read->num_mismatch();
   if(_right_read)
      edits += _right_read->num_mismatch();
   return edits;
}

int PairedHit::numHits() const
{
   int num = 0;
   if(_left_read){
      num += _left_read->numHits();
   }
   if(_right_read){
      num += _right_read->numHits();
   }
   return num;
}

bool PairedHit::is_multi() const
{
   return numHits() > 1;
}

bool PairedHit::contains_splice() const
{
   if(_right_read)
      return (_left_read->contains_splice() || _right_read->contains_splice());
   else
      return _left_read->contains_splice();
}

ReadID PairedHit::read_id() const
{
   if(_left_read) return _left_read->read_id();
   if(_right_read) return _right_read->read_id();
   return 0;
}

RefID PairedHit::ref_id() const
{
   if(_left_read && _right_read)
      assert(_left_read->interval().seq_id() == _right_read->interval().seq_id());
   if(_left_read)
      return _left_read->interval().seq_id();
   if(_right_read)
      return _right_read->interval().seq_id();
   return -1;
}

bool PairedHit::operator==(const PairedHit& rhs)
{
   if((rhs.left_read() == nullptr) != (left_read() == nullptr))
      return false;
   if((rhs.right_read() == nullptr) != (right_read() == nullptr))
      return false;
   if(left_read()){
      if(left_read() != rhs.left_read()) return false;
   }
   if(right_read()){
      if(right_read() != rhs.right_read()) return false;
   }
   return true;
}

RefID RefSeqTable::get_id(const string& name) {
   if (name == "*") return -1;
   unordered_map<string,int>::const_iterator it = _name2id.find(name);
   if(it != _name2id.end()) return it->second;
   else {
      int id = _name2id.size();
      _name2id.insert(make_pair(name, id));
      _id2name.resize(_name2id.size());
      _id2name[id] = name;
      return id;
   }
   return 0;
}

