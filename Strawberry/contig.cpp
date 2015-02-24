/*
 * contig.cpp
 *
 *  Created on: Jan 19, 2015
 *      Author: ruolin
 */



#include "contig.h"
#include <algorithm>
#include <iostream>
using namespace std;

GenomicFeature::GenomicFeature(const Op_t& cc, uint offset, int len):
         _code(cc),
         _genomic_offset(offset),
         _genomic_length(len)
{
   assert (_genomic_length >= 0);
}

void GenomicFeature::left(uint left)
{
      int right = _genomic_offset + _genomic_length;
      _genomic_offset = left;
      _genomic_length = right -left;
}

uint GenomicFeature::left() const { return _genomic_offset;}

uint GenomicFeature::right() const
{
   return _genomic_offset + _genomic_length-1;
}

void GenomicFeature::right(uint right)
{
      _genomic_length = right - _genomic_offset + 1;
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
   if(left_off + op._genomic_length > left && left_off < right){
      if(left_off > left){
         if(left_off + op._genomic_length <= right + 1)
            len += op._genomic_length;
         else
            len += right -left_off;
      }
      else{
         if(left_off + op._genomic_length <= right +1)
            len += (left_off + op._genomic_length - left);
         else
            return right - left;
      }
   }
   return len;
}

bool GenomicFeature::operator==(const GenomicFeature & rhs) const
{
   return ( _code == rhs._code &&
         _genomic_offset == rhs._genomic_offset &&
         _genomic_length == rhs._genomic_length);
}

bool GenomicFeature::operator<(const GenomicFeature & rhs) const
{
   if(_genomic_offset != rhs._genomic_offset)
      return _genomic_offset < rhs._genomic_offset;
   if(_genomic_length != rhs._genomic_length)
      return _genomic_length < rhs._genomic_length;
   return false;
}

Contig::Contig(RefID ref_id, char strand, const vector<GenomicFeature> &feats, bool is_ref):
      _genomic_feats(feats),
      _is_ref(is_ref),
      _ref_id(ref_id),
      _strand(strand)
   {
      assert(_genomic_feats.front()._code == S_MATCH);
      assert(_genomic_feats.back()._code == S_MATCH);
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

bool Contig::operator<(const Contig &rhs) const
{
   if(_ref_id != rhs._ref_id){
      return _ref_id < rhs._ref_id;
   }
   if(left() != rhs.left()){
      return left() < rhs.left();
   }
   return false;
}

RefID Contig::ref_id() const
{
   return _ref_id;
}

const char Contig::strand() const
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

bool Contig::overlaps_directional(const Contig &lhs, const Contig &rhs){
   if(lhs.ref_id() != rhs.ref_id())
      return false;
   if(lhs.strand() != rhs.strand())
      return false;
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
