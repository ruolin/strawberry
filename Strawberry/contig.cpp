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
GenomicInterval::GenomicInterval(string chr, uint l, uint r, char o) {
      _chrom = chr;
      transform(chr.begin(), chr.end(), chr.begin(), ::tolower);
      if (l>r) { _left = r; _right = l;}
      else { _left = l; _right = r;}
      switch(o){
      case '+':
         _strand = kStrandPlus;
         break;
      case '-':
         _strand = kStrandMinus;
         break;
      case '.':
         _strand = kStrandUnknown;
         break;
      default:
         cerr<<"Error: stand must be '+', '-' or '.'"<<endl;
         break;
      }
   }



uint GenomicInterval::left() const{ return _left;}

uint GenomicInterval::right() const { return _right;}

void GenomicInterval::set_left(uint l) {_left = l;}

void GenomicInterval::set_right(uint r) {_right = r;}

char GenomicInterval::strand() const { return _strand;}

string GenomicInterval::chrom() const { return _chrom;}

uint GenomicInterval::len() const { return _right-_left+1;}


bool GenomicInterval::overlap(const GenomicInterval& other, bool nonStrandness) const
{
     if( _chrom != other._chrom) return false;
     if( !nonStrandness && other._strand != kStrandUnknown && _strand != kStrandUnknown && other._strand != _strand) return false;
     return _left < other._left ? ( other._left <= _right) : (_left <= other._right);
}

bool GenomicInterval::isContainedIn(const GenomicInterval &other, bool nonStrandness) const
{
     if( other._chrom != _chrom) return false;
     if( !nonStrandness && other._strand != kStrandUnknown && _strand != kStrandUnknown && other._strand != _strand) return false;
     if (_left < other._left || _right > other._right) return false;
     return true;
}

bool GenomicInterval::contain(const GenomicInterval &d, bool nonStrandness) const
{
     return d.isContainedIn(*this, nonStrandness);
}

  //return the length of overlap between two segments
uint GenomicInterval::overlapLen(const GenomicInterval& other) const
{
     if (!other.overlap(*this)) cerr<<"this two interval does not overlap"<<endl;
     if (_left<other._left) {
        if (other._left>_right) return 0;
        return (other._right>_right) ? _right-other._left+1 : other._right-other._left+1;
        }
       else { //r->start<=start
        if (_left>other._right) return 0;
        return (other._right<_right)? other._right-_left+1 : _right-_left+1;
        }
}

bool GenomicInterval::operator==(const GenomicInterval& rhs) const
{
     if ( rhs._chrom != _chrom) return false;
     if( rhs._strand != kStrandUnknown && _strand != kStrandUnknown && rhs._strand != _strand) return false;
     return (_left == rhs._left && _right == rhs._right);
}

bool GenomicInterval::operator>(const GenomicInterval& rhs) const
{
     if ( rhs._chrom != _chrom) cerr<<"cannot compare for different chrom"<<endl;
     return (_left==rhs._left)?(_right>rhs._right):(_left>rhs._left);
}

bool GenomicInterval::operator<(const GenomicInterval& rhs) const
{
     if ( rhs._chrom != _chrom) cerr<<"cannot compare for different chrom"<<endl;
     return (_left == rhs._left)?(_right < rhs._right):(_left < rhs._left);
}

GenomicFeature::GenomicFeature(FeatType ftype, GenomicInterval iv):
      _feature(ftype), _iv(iv){};

