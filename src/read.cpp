/*
>HEADER
    Copyright (c) 2015 Ruolin Liu rliu0606@gmail.com
    This file is part of Strawberry.
    Strawberry is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Strawberry is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Strawberry.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
*/


#include <algorithm>
#include<assert.h>
#include <stdexcept>
#include <string.h>
//#include <typeinfo>
#include<cxxabi.h>
#include "read.hpp"
//#include "kmer.h"

using namespace std;
void mean_and_sd_insert_size(const vector<int> & vec, double & mean, double &sd){
   double sum = accumulate(vec.begin(), vec.end(), 0.0);
   mean = sum / vec.size();
   double sq_sum = inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
   sd = std::sqrt(sq_sum / vec.size() - mean * mean);

}


ReadHit::ReadHit(
   ReadID readID,
   GenomicInterval iv,
   const vector<CigarOp> & cigar,
   RefID partnerRef,
   int partnerPos,
   int numMismatch,
   int numHit,
   uint32_t samFlag,
   double mass,
   char* seq):
      _read_id(readID),
      _iv(iv),
      _cigar(cigar),
      _partner_ref_id(partnerRef),
      _partner_pos(partnerPos),
      _num_mismatch(numMismatch),
      _num_hits(numHit),
      _sam_flag(samFlag)
{
   if(seq != NULL){
      _seq = string(seq);
   }

   if(is_singleton()){
      _read_mass = 1.0/_num_hits;
   }else{
      _read_mass =  0.5/_num_hits;
   }
}

const vector<CigarOp>& ReadHit::cigar() const
{
   return _cigar;
}

uint ReadHit::read_len() const
{
   uint len = 0;
   for(size_t i =0; i< _cigar.size(); ++i){
      switch(_cigar[i]._type)
      {
      case MATCH:
      case SOFT_CLIP:
      case INS:
         len +=_cigar[i]._length;
         break;
      default:
         break;
      }
   }
   return len;
}

uint ReadHit::intron_len() const
{
   uint len = 0;
   for(size_t i=0; i< _cigar.size(); ++i){
      if(_cigar[i]._type == REF_SKIP)
         len += _cigar[i]._length;
   }
   return len;
}

vector<pair<uint,uint>> ReadHit::intron_coords() const
  /*
    *  Only be called if objects contains intron
    */
{
   assert(contains_splice());
   vector<pair<uint,uint>> coords;
   uint start = left();
   uint end = 0;
   for(size_t i=0; i< _cigar.size(); ++i){
      if(_cigar[i]._type != REF_SKIP)
         start += _cigar[i]._length;
      else{
         end = start + _cigar[i]._length -1;
         coords.emplace_back(start,end);
         start += _cigar[i]._length;
      }
   }
   return coords;
}

double ReadHit::mass() const
{
   return _read_mass;
}

double ReadHit::raw_mass() const
{
   if(is_singleton()){
      return  1.0/_num_hits;
   }else{
      return  0.5/_num_hits;
   }
}

void ReadHit::mass(double m)
{
   _read_mass = m;
}

bool ReadHit::contains_splice()const{
      for (size_t i = 0; i < _cigar.size(); ++i){

            if (_cigar[i]._type == REF_SKIP)
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

Strand_t ReadHit::strand() const {return _iv.strand();}

//vector<CigarOp> ReadHit::cigars() const {
//   return _cigar;
//}

uint ReadHit::left() const { return _iv.left();}

int ReadHit::numHits() const { return _num_hits;}
int ReadHit::num_mismatch() const { return _num_mismatch;}

bool ReadHit::is_singleton() const
{
//#ifdef DEBUG
//   cout<<partner_pos()<<endl;
//#endif
   return (partner_pos() == 0 ||
         partner_ref_id() == -1 ||
         partner_ref_id() != ref_id());
}

bool ReadHit::reverseCompl() const
{
   return _sam_flag  & 0x10;
}

// not considering read orientation and cigar string
bool ReadHit::operator<(const ReadHit& rhs) const
{

   if(left() == rhs.left())
      return right() < rhs.right();
   else
      return left() < rhs.left();
}

bool ReadHit::operator==(const ReadHit& rhs) const
{
   //assert(!_cigar.empty() && !rhs._cigar.empty());
   //return (_iv == rhs.interval());
   if(left() == rhs.left())
      return right() == rhs.right();
   else
      return false;
}

bool ReadHit::operator!=(const ReadHit& rhs) const
{
   return !(*this == rhs);
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

InsertSize::InsertSize():
                     _mean(kInsertSizeMean),
                     _sd(kInsertSizeSD),
                     _use_emp(false){}

InsertSize::InsertSize(double mean, double sd):
                     _mean(mean),
                     _sd(sd),
                     _use_emp(false){}

InsertSize::InsertSize(const vector<int> frag_lens):_use_emp(true)
{
   _total_reads = frag_lens.size();
   mean_and_sd_insert_size(frag_lens, _mean, _sd);
   auto result = minmax_element(frag_lens.begin(), frag_lens.end());
   if (verbose){
      cerr<<"Calculated averaged fragment length is: "<<_mean<<endl;
      cerr<<"Calculated fragment length sd is: "<<_sd<<endl;
      cerr<<"Min fragment length is: "<<*result.first<<endl;
      cerr<<"Max fragment length is: "<<*result.second<<endl;
   }
   _start_offset = *result.first;
   _end_offset = *result.second;
   _emp_dist.resize(_end_offset - _start_offset +1, 0);
   for(size_t i = 0; i< frag_lens.size(); ++i){
      _emp_dist[frag_lens[i]-_start_offset] ++ ;
   }

#ifdef DEBUG
   //cout<<"number of fragments: "<<frag_lens.size()<<endl;
   size_t n = accumulate(_emp_dist.begin(), _emp_dist.end(), 0);
   assert(n == _total_reads);
   //cout<<"empirical distribution len: "<<n<<endl;
#endif
}

bool InsertSize::empty() const
{
   return _emp_dist.empty();
}

double InsertSize::emp_dist_pdf(uint insert_size) const
{
   if(_use_emp){
      double ret = 0.0;
      if (insert_size < _start_offset || insert_size > _end_offset) {
      } else {
         ret = _emp_dist[insert_size - _start_offset] / _total_reads;
      }

      if (ret == 0.0) {
         double p = normal_pdf( (double) insert_size, _mean, _sd);
         if(p > 0) return p;
         else return 0.0;
      } else {
         return ret;
      }
   }

   else{
      double p = normal_pdf( (double) insert_size, _mean, _sd);
      if(p > 0) return p;
      else return 0.0;
   }
}

//double InsertSize::truncated_normal_pdf(uint insert_size) const
//{
//   using boost::math::normal;
//   normal standard_normal;
//   double numerator = 1/_sd * pdf(standard_normal, (insert_size - _mean)/_sd);
//   double denominator = 1 - cdf(standard_normal, (0 - _mean)/_sd);
//   assert(denominator != 0);
//   return numerator/denominator;
//}


HitFactory::HitFactory(ReadTable &reads_table, RefSeqTable &ref_table, string hit_file_name):
   _reads_table(reads_table),_ref_table(ref_table), _hit_file_name(hit_file_name){}

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
            //str2lower(fields[1]);
            RefID _ID = _ref_table.get_id(fields[1]);
            if(_ID != _num_seq_header_recs-1){
               cerr<<"Sort order of reads in BAM not consistent."<<endl;
               exit(0);
            }
         }
      }
   }

   if(cols[0] == "@RG"){
      for(auto &i: cols){
         vector<string> fields;
         split(i, ":", fields);
         if(fields[0] == "PL"){
            platform_t p = str2platform(fields[1]);
            assert(_assay_props._platform == UNKNOWN_PLATFORM);
            _assay_props._platform = p;
         }
      }
   }
   return true;
}


BAMHitFactory::BAMHitFactory(const string& bam_file_name,
                            ReadTable &read_table,
                            RefSeqTable &ref_table):
                 HitFactory(read_table, ref_table, bam_file_name)
{
   _hit_file = samopen(bam_file_name.c_str(), "rb", 0);
   memset(&_next_hit, 0, sizeof(_next_hit));
   if(_hit_file == NULL || _hit_file->header == NULL){
      cerr<<"Fail to open BAM file..."<<endl;
      exit(0);
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
      LOG_WARN("No BAM header");
      return false;
   }

   if(header->l_text >= MAX_HEADER_LEN ){
      LOG_WARN("BAM header too large");
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

   /*For unmapped read or partner
    * core.pos = -1 or
    * core.mpos = -1
    * */
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
      //str2lower(mate_text_name);
   }
   RefID parterner_ref_id = ref_table().get_id(mate_text_name);

   ReadID readid = HitFactory::reads_table().get_id(bam1_qname(hit_buf));
   vector<CigarOp> cigar;
   bool is_spliced_alignment = false;
   int num_hits = 1;
   if( (sam_flag & 0x4) || target_id < 0 ){ // unmapped reads
         return false;
      bh = ReadHit(readid,
                   GenomicInterval(),
                   cigar,
                   parterner_ref_id,
                   0,
                   0,
                   num_hits,
                   sam_flag,
                   0.0,
                   NULL
                   );
      return true;
   }

   string text_name = _hit_file->header->target_name[target_id];
   //str2lower(text_name);
   RefID ref_id = ref_table().get_id(text_name);

   if(target_id >= _hit_file->header->n_targets){
      fprintf(stderr, "BAM error: file contains hits to sequences not in BAM file header");
      return false;
   }

   int read_len = 0;
   int eff_read_len = 0;
   for (int i=0; i<hit_buf->core.n_cigar; ++i){
      int length = bam1_cigar(hit_buf)[i] >> BAM_CIGAR_SHIFT;
      if(length <= 0){
         fprintf(stderr, "BAM error: CIGAR op has zero length (%s)\n",bam1_qname(hit_buf));
         return false;
      }
      CigarOpCode _type;
      switch(bam1_cigar(hit_buf)[i] & BAM_CIGAR_MASK){
      case BAM_CMATCH: _type = MATCH;
         read_len += length;
         eff_read_len += length;
         cigar.push_back(CigarOp(_type, length));
         break;
      case BAM_CINS: _type = INS; // INSERTION does not increase read length
         cigar.push_back(CigarOp(_type, length));
         break;
      case BAM_CDEL: _type = DEL;
         read_len += length;
         cigar.push_back(CigarOp(_type, length));
         break;
      case BAM_CSOFT_CLIP:
         _type = SOFT_CLIP;
         cigar.push_back(CigarOp(_type, length));
         break;
      case BAM_CHARD_CLIP:
         _type = HARD_CLIP;
         break;
      case BAM_CPAD:
         _type = PAD;
         break;
      case BAM_CREF_SKIP:
         _type = REF_SKIP;
         is_spliced_alignment = true;
         read_len += length;
         cigar.push_back(CigarOp(_type, length));
         if(length > kMaxIntronLength || length < kMinIntronLength){
            LOG_ERR("At read ", bam1_qname(hit_buf), " has unreasonable intron size: ", length);
            return false;
         }
         break;
      default: return false;
      }
   }

   /*Filtering based on Cigar.
    * DEL and INS must be sandwiched by MATCH
    * */
   for(int i=0; i != cigar.size(); ++i){
      if(cigar[i]._type == INS || cigar[i]._type == DEL){
         if(i-1 <= 0 || i+1 >= cigar.size())
            return false;
         if(cigar[i-1]._type != MATCH || cigar[i+1]._type != MATCH)
            return false;
      }
   }

   if(eff_read_len <= 1 ) return false;

   string mate_ref;

   if(sam_flag & BAM_FPAIRED) //paired-end read
   {
      SINGLE_END_EXP = false;
      if(mate_target_id == target_id){
         mate_ref = _hit_file->header->target_name[mate_target_id];
      }
      else {
         if(sam_flag & 0x8){ // next segment unmapped
            if( mate_pos != 0 || mate_target_id >= 0) {
               cerr<<"read "<<bam1_qname(hit_buf)<<" is ill-formed"<<endl;
            }
         }
         else{
            LOG_WARN("At read ", bam1_qname(hit_buf), " read pair aligns to different chromosome");
            return false;
         }

      }
   }
   else{ // single-end read
      if( mate_pos != 0 || mate_target_id >= 0) {
         cerr<<"read "<<bam1_qname(hit_buf)<<" is ill-formed"<<endl;
      }
   }

   Strand_t source_strand = Strand_t::StrandUnknown;
   unsigned char num_mismatches = 0;

   uint8_t* ptr = bam_aux_get(hit_buf, "XS");
   if(ptr){
      char src_strand_char = bam_aux2A(ptr);
      switch(src_strand_char){
      case '+':
         source_strand = Strand_t::StrandPlus;
         break;
      case '-':
         source_strand = Strand_t::StrandMinus;
         break;
      default:
         //LOG_WARN("At read ", bam1_qname(hit_buf), " parsing spliced alignment without known strand information");
         break;
      }
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

//   if(is_spliced_alignment){
//      if(source_strand == Strand_t::StrandUnknown)
//         fprintf(stderr, "BAM record error: Unknown strand for spliced alignment, XS attribute is missing\n");
//   }

   if(use_only_unique_hits && num_hits > 1)
      return false;
//   if(use_only_paired_hits && ( sam_flag & 0x8 || mate_target_id != target_id ))
//      return false;

   /*
    * Get read sequencing
    * */
   char *qseq = (char *) malloc(hit_buf->core.l_qseq+1);
   unsigned char *s   = bam1_seq(hit_buf);
   int n= 0;
   for(n = 0;n<(hit_buf->core.l_qseq);n++) {
     int v = bam1_seqi(s,n);
     qseq[n] = bam_nt16_rev_table[v];
   }
   qseq[n] = 0;
   bh = ReadHit(
               readid,
               GenomicInterval(ref_id, pos, pos+read_len-1, source_strand),
               cigar,
               parterner_ref_id,
               mate_pos,
               num_mismatches,
               num_hits,
               sam_flag,
               mass,
               qseq
               );
   free(qseq);
   qseq = NULL;
   return true;
}


PairedHit::PairedHit(ReadHitPtr leftRead, ReadHitPtr rightRead):
            _left_read(move(leftRead)), _right_read(move(rightRead))
{
}


void PairedHit::weighted_mass(double m)
{
   _mass = m;
}

double PairedHit::weighted_mass() const
{
   return _mass;
}

void PairedHit::init_raw_mass()
{
   assert(_mass == 0.0);
   if(_left_read)
      _mass += _left_read->mass();
   if(_right_read)
      _mass += _right_read->mass();
}

//void PairedHit::set_kmers(int num_kmers )
//{
//   if(_left_read ){
//      for(int i = 0; i< num_kmers; ++i){
//         string sub = _left_read->_seq.substr(i, Kmer::_k);
//         Kmer k(sub.c_str());
//         _left_kmers.push_back(k);
//      }
//   }

//   if(_right_read){
//      string rev_seq = _right_read->_seq;
//      reverse(rev_seq.begin(), rev_seq.end());
//      for(int i = 0; i< num_kmers; ++i){
//         string sub = rev_seq.substr(i, Kmer::_k);
//         Kmer k(sub.c_str());
//         _right_kmers.push_back(k);
//      }
//   }

//}

const ReadHit& PairedHit::left_read_obj() const {return *_left_read;}

Strand_t PairedHit::strand() const {
   if(_right_read && _left_read){
      assert(_right_read->strand() == _left_read->strand() ||
            _right_read->strand() == Strand_t::StrandUnknown ||
            _left_read->strand() == Strand_t::StrandUnknown);
      if(_left_read->strand() != Strand_t::StrandUnknown)
         return _left_read->strand();
      else
         return _right_read->strand();
   } else if (_left_read)
      return _left_read->strand();
   else if (_right_read)
      return _right_read->strand();
   else
      assert(false);
}

void PairedHit::set_left_read(ReadHitPtr lr)
{
   _left_read = move(lr);
}

const ReadHit& PairedHit::right_read_obj() const {return *_right_read;}

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
   bool l_res, r_res;
   if(_left_read)
      l_res = _left_read->contains_splice();
   else
      l_res = false;
   if(_right_read)
      r_res = _right_read->contains_splice();
   else
      r_res = false;

   return l_res || r_res;
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

double PairedHit::raw_mass() const
{
   double m = 0.0;
   if(_left_read)
      m += _left_read->mass();
   if(_right_read)
      m += _right_read->mass();
   return m;
}

bool PairedHit::operator==(const PairedHit& rhs) const
{
   if((rhs._left_read == nullptr) != (_left_read == nullptr))
      return false;
   if((rhs._right_read == nullptr) != (_right_read == nullptr))
      return false;
   if(_left_read){
      if(left_read_obj() != rhs.left_read_obj()) return false;
   }
   if(_right_read){
      if(right_read_obj() != rhs.right_read_obj()) return false;
   }
   return true;
}

bool PairedHit::operator !=(const PairedHit& rhs) const
{
   return !(*this == rhs);
}

bool PairedHit::operator<(const PairedHit& rhs) const
{
   if(left_pos() == rhs.left_pos())
      return right_pos() < rhs.right_pos();
   else
      return left_pos() < rhs.left_pos();
}

//bool PairedHit::operator<(const PairedHit& rhs) const
//{
//   if(_left_read != nullptr && rhs._left_read != nullptr){
//      if(left_read_obj() == rhs.left_read_obj()){
//         if(_right_read ==nullptr && rhs._right_read != nullptr) return true;
//         else if(_right_read != nullptr && rhs._right_read ==nullptr) return false;
//         else if(_right_read == nullptr && rhs._right_read == nullptr) return false;
//         else return right_read_obj() < rhs.right_read_obj();
//      }
//      else{
//         return left_read_obj() < rhs.left_read_obj();
//      }
//   }
//   else if(_left_read == nullptr && rhs._left_read ==nullptr){
//      return right_read_obj() < rhs.right_read_obj();
//   }
//   else if(_left_read == nullptr && rhs._left_read !=nullptr){
//      if(rhs._right_read != nullptr){
//         return right_read_obj() < rhs.right_read_obj();
//      }
//      else{
//         if(right_read_obj() == rhs.left_read_obj()) return true;
//         return right_read_obj() < rhs.left_read_obj();
//      }
//   }
//   else{ // _left_read != nullptr && rhs._left_read ==nullptr
//      if(_right_read != nullptr)
//         return right_read_obj() < rhs.right_read_obj();
//      else
//         return left_read_obj() < rhs.right_read_obj();
//   }
//}

RefID RefSeqTable::get_id(string& name) {
   if (name == "*") return -1;
   string raw_name = name;
   str2lower(name);
   unordered_map<string,int>::const_iterator it = _name2id.find(name);
   if(it != _name2id.end()) return it->second;
   else {
      int id = _name2id.size();
      _name2id.insert(make_pair(name, id));
      _id2name.resize(_name2id.size());
      _id2name[id] = name;
      _id_2_real_name.resize(_name2id.size());
      _id_2_real_name[id] = raw_name;
      return id;
   }
   return 0;
}

const string RefSeqTable::ref_real_name(int id) const{
   return _id_2_real_name[id];
}

void PairedHit::add_2_collapse_mass(double add){
   _collapse_mass += add;
}

double PairedHit::collapse_mass() const {
   return _collapse_mass;
}
