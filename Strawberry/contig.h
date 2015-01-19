/*
 * contig.h
 *
 *  Created on: Jan 19, 2015
 *      Author: ruolin
 */

#ifndef CONTIG_H_
#define CONTIG_H_

#include<string>
#include "GBase.h"

static const char strand_plus = '+';
static const char strand_minus = '-';
static const char strand_unknown = '.';

class GenomicInterval {
private:
   uint _left; //start<end always!
   uint _right;
   std::string _chrom;
   char _strand;


public:

   GenomicInterval()=default;
   GenomicInterval(const std::string chr, uint l, uint r, char o);

  uint left() const;
  uint right() const;
  void set_left(uint l);
  void set_right(uint r);
  const char strand() const;
  const std::string chrom() const;
  uint len() const;

  //check for overlap with other segment
  bool overlap(const GenomicInterval &other, bool nonStrandness = true) const;

  bool isContainedIn(const GenomicInterval &other, bool nonStrandness = true) const;

  bool contain(const GenomicInterval &other, bool nonStrandness = true) const;

  //return the length of overlap between two segments
  uint overlapLen(const GenomicInterval& other) const;
  //comparison operators required for sorting
  bool operator==(const GenomicInterval& d) const;
  bool operator>(const GenomicInterval& d) const;
  bool operator<(const GenomicInterval& d) const;
};

class GenomicFeature{
// Represent either exon or intron or unknown gap.
public:
   enum FeatType{FEAT_EXON, FEAT_INTRON, FEAT_GAP};
   GenomicInterval _iv;
   FeatType _feature;
public:
   GenomicFeature() = default;
   GenomicFeature(FeatType ftype, GenomicInterval iv);
   GenomicFeature(FeatType ftype, const std::string chr, uint l, uint r, char o);
};


#endif /* CONTIG_H_ */
