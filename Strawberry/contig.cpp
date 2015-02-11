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

void GenomicFeature::g_left(uint left)
{
      int right = _genomic_offset + _genomic_length;
      _genomic_offset = left;
      _genomic_length = right -left;
}

uint GenomicFeature::g_left() const { return _genomic_offset;}

uint GenomicFeature::g_right() const
{
   return _genomic_offset + _genomic_length-1;
}

void GenomicFeature::g_right(uint right)
{
      _genomic_length = right - _genomic_offset + 1;
}

bool GenomicFeature::overlap_in_genome(const GenomicFeature& lhs, const GenomicFeature& rhs)
{
   if (lhs.g_left() >= rhs.g_left() && lhs.g_left() < rhs.g_right())
      return true;
   if (lhs.g_right() > rhs.g_left() && lhs.g_right() < rhs.g_right())
      return true;
   if (rhs.g_left() >= lhs.g_left() && rhs.g_left() < lhs.g_right())
      return true;
   if (rhs.g_right() > lhs.g_left() && rhs.g_right() < lhs.g_right())
      return true;
   return false;
}

bool GenomicFeature::contains(const GenomicFeature& other) const
{
   if (g_left() <= other.g_left() && g_right() >= other.g_right())
      return true;
   return false;
}

bool GenomicFeature::properly_contains(const GenomicFeature& other) const
{
   if( (g_left() < other.g_left() && g_right() >= other.g_right() ) ||
         (g_left() <= other.g_left() && g_right() > other.g_right()) )
      return true;
   return false;
}

int match_length(const GenomicFeature &op, int left, int right)
{
   int len = 0;
   int left_off = op.g_left();
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

Contig::Contig(RefID ref_id, char strand, vector<GenomicFeature> &feats, bool is_ref):
      _genomic_feats(move(feats)),
      _is_ref(is_ref),
      _ref_id(ref_id),
      _strand(strand)
   {
      assert(_genomic_feats.front()._code == S_MATCH);
      assert(_genomic_feats.back()._code == S_MATCH);
   }

