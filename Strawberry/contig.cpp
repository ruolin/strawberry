/*
 * contig.cpp
 *
 *  Created on: Jan 19, 2015
 *      Author: ruolin
 */



#include "contig.h"
#include <algorithm>
#include <iostream>
#include "read.h"
using namespace std;



bool readhit_2_genomicFeats(const ReadHit & rh, vector<GenomicFeature> & feats){
   vector<CigarOp> cig = rh.cigar();
   uint offset = rh.left();
   assert(!cig.empty());
   for(size_t i = 0; i< cig.size(); ++i){
      switch(cig[i]._type)
      {
      case MATCH:
         feats.push_back(GenomicFeature(Match_t::S_MATCH, offset, cig[i]._length));
         offset += cig[i]._length;
         break;
      case REF_SKIP:
         feats.push_back(GenomicFeature(Match_t::S_INTRON, offset, cig[i]._length));
         offset += cig[i]._length;
         break;
      case DEL:
         if(i<1 || i+1 == cig.size() || cig[i-1]._type != MATCH || cig[i+1]._type != MATCH){
            LOG_INPUT("Read at reference id: ", rh.ref_id()+1, " and position ", rh.left(), " has suspicious DELETION");
            return false;
         }
         feats.back()._match_op._len += cig[i]._length;
         offset += cig[i]._length;
         ++i;
         feats.back()._match_op._len += cig[i]._length;
         offset += cig[i]._length;
         break;
      case INS:
         if(i<1 || i+1 == cig.size() || cig[i-1]._type != MATCH || cig[i+1]._type != MATCH){
            LOG_INPUT("Read at reference id: ", rh.ref_id()+1, " and position ", rh.left(), " has suspicious INSERTION");
            return false;
         }
         ++i;
         feats.back()._match_op._len += cig[i]._length;
         offset += cig[i]._length;
         break;
      case SOFT_CLIP:
         break;
      default:
         LOG_INPUT("Read at reference id: ", rh.ref_id()+1, " and position ", rh.left(), " has unknown cigar");
         return false;
      }
   }
   return true;
}

GenomicFeature::GenomicFeature(const Match_t& code, uint offset, int len):
         _genomic_offset(offset),
         _match_op (MatchOp(code, (uint32_t)len))
{
   assert (len >= 0);
}

void GenomicFeature::left(uint left)
{
      int right = _genomic_offset + _match_op._len;
      _genomic_offset = left;
      _match_op._len = right -left;
}

uint GenomicFeature::left() const { return _genomic_offset;}

uint GenomicFeature::right() const
{
   return _genomic_offset + _match_op._len-1;
}

void GenomicFeature::right(uint right)
{
      _match_op._len = right - _genomic_offset + 1;
}

bool GenomicFeature::overlap_in_genome(const GenomicFeature& lhs, const GenomicFeature& rhs)
{
   if (lhs.left() >= rhs.left() && lhs.left() < rhs.right())
      return true;
   if (lhs.right() > rhs.left() && lhs.right() < rhs.right())
      return true;
   if (rhs.left() >= lhs.left() && rhs.left() < lhs.right())
      return true;
   if (rhs.right() > lhs.left() && rhs.right() < lhs.right())
      return true;
   return false;
}

bool GenomicFeature::contains(const GenomicFeature& other) const
{
   if (left() <= other.left() && right() >= other.right())
      return true;
   return false;
}

bool GenomicFeature::properly_contains(const GenomicFeature& other) const
{
   if( (left() < other.left() && right() >= other.right() ) ||
         (left() <= other.left() && right() > other.right()) )
      return true;
   return false;
}

int match_length(const GenomicFeature &op, int left, int right)
{
   int len = 0;
   int left_off = op.left();
   if(left_off + op._match_op._len > left && left_off < right){
      if(left_off > left){
         if(left_off + op._match_op._len <= right + 1)
            len += op._match_op._len;
         else
            len += right -left_off;
      }
      else{
         if(left_off + op._match_op._len <= right +1)
            len += (left_off + op._match_op._len - left);
         else
            return right - left;
      }
   }
   return len;
}

bool GenomicFeature::operator==(const GenomicFeature & rhs) const
{
   return ( _match_op._code == rhs._match_op._code &&
         _genomic_offset == rhs._genomic_offset &&
         _match_op._len == rhs._match_op._len);
}

bool GenomicFeature::operator<(const GenomicFeature & rhs) const
{
   if(_genomic_offset != rhs._genomic_offset)
      return _genomic_offset < rhs._genomic_offset;
   if(_match_op._len != rhs._match_op._len)
      return _match_op._len < rhs._match_op._len;
   return false;
}

Contig::Contig(RefID ref_id, Strand_t strand, const vector<GenomicFeature> &feats, double mass, bool is_ref):
      _genomic_feats(feats),
      _is_ref(is_ref),
      _ref_id(ref_id),
      _strand(strand),
      _mass(mass)
   {
      assert(_genomic_feats.front()._match_op._code == Match_t::S_MATCH);
      assert(_genomic_feats.back()._match_op._code == Match_t::S_MATCH);
   }

Contig::Contig(const PairedHit& ph):
      _is_ref(false),
      _ref_id(ph.ref_id()),
      _strand(ph.strand())
{
   vector<GenomicFeature> g_feats;
   if(ph._left_read){
      if(readhit_2_genomicFeats(ph.left_read_obj(), g_feats)){
         if(ph._right_read){
            int gap_len = (int)ph._right_read->left() - (int)ph._left_read->right() -1 ;
            if( gap_len < 0){
               gap_len = 0;
               LOG_INPUT("Read at reference id: ", ph.ref_id()+1, " and position ", ph.left_pos(), " has suspicious GAP");
            }
            g_feats.push_back(GenomicFeature(Match_t::S_GAP, ph._left_read->right()+1, (uint)gap_len));
            readhit_2_genomicFeats(ph.right_read_obj(), g_feats);
         }
      }

   }

   else{
      if(ph._right_read){
         readhit_2_genomicFeats(ph.right_read_obj(), g_feats);
      }
      assert(false);
   }

   if(!is_sorted(g_feats.begin(), g_feats.end()))
   {
      sort(g_feats.begin(), g_feats.end());
   }

   uint check_right = ph.left_pos();
   for(auto i = g_feats.cbegin(); i != g_feats.cend(); ++i){
      check_right += i->_match_op._len;
   }
   assert(check_right = ph.right_pos()+1);
   _genomic_feats = move(g_feats);
   _mass = ph.collapse_mass();
}

uint Contig::left() const
{
   //cout<<_genomic_feats.front().left()<<endl;

   return _genomic_feats.front().left();
}

uint Contig::right() const
{

   return _genomic_feats.back().right();
}

float Contig::mass() const
{
   return _mass;
}

void Contig::mass(float m)
{
   _mass = m;
}

bool Contig::operator<(const Contig &rhs) const
{
   if(_ref_id != rhs._ref_id){
      return _ref_id < rhs._ref_id;
   }
   if(left() != rhs.left()){
      return left() < rhs.left();
   }
   if(right() != rhs.right()){
      return right() < rhs.right();
   }
   if(_genomic_feats.size() != rhs._genomic_feats.size()){
      return _genomic_feats.size() < rhs._genomic_feats.size();
   }
   for(size_t i=0; i<_genomic_feats.size(); ++i){
      if(_genomic_feats[i] != rhs._genomic_feats[i])
         return _genomic_feats[i] < rhs._genomic_feats[i];
   }
   return false;
}

RefID Contig::ref_id() const
{
   return _ref_id;
}

Strand_t Contig::strand() const
{
   return _strand;
}

const string Contig::annotated_trans_id() const{
   return _annotated_trans_id;
}

void Contig::annotated_trans_id(string str){
   assert(!str.empty());
   _annotated_trans_id = str;
}

size_t Contig::featSize() const{
   return _genomic_feats.size();
}

bool Contig::overlaps_directional(const Contig &lhs, const Contig &rhs){
   if(lhs.ref_id() != rhs.ref_id())
      return false;
   if(lhs.strand() != rhs.strand())
      return false;
   return overlaps_locally(lhs.left(), lhs.right(), rhs.left(), rhs.right());
}


