/*
 * common.cpp
 *
 *  Created on: Feb 2, 2015
 *      Author: ruolin
 */
#include "common.h"
#include <algorithm>

static char msg[4096];
int fileExists(const char* fname)
{
  struct stat stFileInfo;
  int r=0;
  // Attempt to get the file attributes
  int fs = stat(fname,&stFileInfo);
  if (fs == 0) {
      r=3;
      // We were able to get the file attributes
      // so the file obviously exists.
      if (S_ISREG (stFileInfo.st_mode)) {
         r=2;
         }
      if (S_ISDIR (stFileInfo.st_mode)) {
          r=1;
          }
      }
  return r;
}

bool endsWith (std::string const &fullString, std::string const &ending)
{
   if (fullString.length() >= ending.length()) {
      return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
   } else {
      return false;
   }
}

bool endsWith(const char* s, const char* suffix)
{
   if (suffix==NULL || s==NULL) return false;
   if (suffix[0]==0) return true; //special case: empty suffix
   int j=strlen(suffix)-1;
   int i=strlen(s)-1;
   if (i<j) return false;
   while (j>=0 && s[i]==suffix[j]) { i--; j--; }
   return (j==-1);
}


void split(const std::string& s, const std::string& delims, std::vector<std::string>& result) \
{
   std::string::size_type lastPos = s.find_first_not_of(delims, 0);
   std::string::size_type pos = s.find_first_of(delims, lastPos);
   while (std::string::npos != pos || std::string::npos != lastPos) {
      result.push_back(s.substr(lastPos, pos - lastPos));
        lastPos = s.find_first_not_of(delims, pos);
        pos = s.find_first_of(delims, lastPos);
    }
}

int stricmp(const char* a, const char* b, int n) {
 if (a==NULL || b==NULL) return a==NULL ? -1 : 1;
 register int ua, ub;
 if (n<0) {
   while ((*a!=0) && (*b!=0)) {
    ua=tolower((unsigned char)*a);
    ub=tolower((unsigned char)*b);
    a++;b++;
    if (ua!=ub) return ua < ub ? -1 : 1;
    }
    return (*a == 0) ? ( (*b == 0) ? 0 : -1 ) : 1 ;
  }
 else {
   while (n && (*a!=0) && (*b!=0)) {
    ua=tolower((unsigned char)*a);
    ub=tolower((unsigned char)*b);
    a++;b++;n--;
    if (ua!=ub) return ua < ub ? -1 : 1;
    }
    //return (*a == 0) ? ( (*b == 0) ? 0 : -1 ) : 1 ;
   if (n==0) return 0;
   else { return (*a == 0) ? ( (*b == 0) ? 0 : -1 ) : 1 ; }
  }
}

void SError(const char* format,...){
  #ifdef __WIN32__
    va_list arguments;
    va_start(arguments,format);
    vsprintf(msg,format,arguments);
    va_end(arguments);
    OutputDebugString(msg);
    fprintf(stderr,"%s",msg); // if a console is available
    MessageBox(NULL,msg,NULL,MB_OK|MB_ICONEXCLAMATION|MB_APPLMODAL);
  #else
    va_list arguments;
    va_start(arguments,format);
    vfprintf(stderr,format,arguments);
    va_end(arguments);
    #ifdef DEBUG
     // modify here if you want a core dump
     abort();
    #endif
  #endif
    exit(1);
}

void SMessage(const char* format,...){
   va_list arguments;
   va_start(arguments,format);
   vsprintf(msg,format,arguments);
   va_end(arguments);
#ifdef __WIN32__
   OutputDebugString(msg);
#endif
   fprintf(stderr,"%s",msg);fflush(stderr);
}

bool SRealloc(pointer* ptr,unsigned long size)
{
   if (size==0) {
      SFree(ptr);
      return true;
   }
   if (*ptr==NULL) {//simple malloc
      void *p=malloc(size);
      if (p != NULL) {
         *ptr=p;
         return true;
      }
      else return false;
   }

   else {//realloc
      void *p=realloc(*ptr,size);
      if (p) {
         *ptr=p;
         return true;
      }
      return false;
   }
}

char* SlineReader::getLine(FILE* stream, off_t& f_pos)
{
   if (pushed) { pushed=false; return buf; }
   //reads a char at a time until \n and/or \r are encountered
   len=0;
   int c=0;
   while ((c=getc(stream))!=EOF) {
      if (len>=allocated-1) {
        allocated+=1024;
        SREALLOC(buf, allocated);
      }
      if (c=='\n' || c=='\r') {
         buf[len]='\0';
         if (c=='\r') { //DOS file -- special case
            if ((c=getc(stream))!='\n') ungetc(c,stream);
            else f_pos++;
         }
         ++f_pos;
         ++lcount;
         return buf;
      }
      f_pos++;
      buf[len++]=(char)c;
   } //while end

   if (c==EOF) {
      isEOF=true;
      if (len==0) return nullptr; // empty file
      buf[len]='\0';
      ++lcount;
      return buf;
   }
}

GenomicInterval::GenomicInterval(int chr, uint l, uint r, Strand_t o) :
   _seq_id(chr),
   _strand(o)
{
      if (l>r) { _left = r; _right = l;}
      else { _left = l; _right = r;}
}



uint GenomicInterval::left() const{ return _left;}

uint GenomicInterval::right() const { return _right;}

void GenomicInterval::set_left(uint l) {_left = l;}

void GenomicInterval::set_right(uint r) {_right = r;}

Strand_t GenomicInterval::strand() const { return _strand;}

int GenomicInterval::seq_id() const { return _seq_id;}

void GenomicInterval::set_seq_id(int id) {
   _seq_id = id;
}
uint GenomicInterval::len() const { return _right-_left+1;}


bool GenomicInterval::overlap(const GenomicInterval& other, bool nonStrandness) const
{
     if( _seq_id != other._seq_id) return false;
     if( !nonStrandness && other._strand != Strand_t::StrandUnknown && _strand != Strand_t::StrandUnknown && other._strand != _strand) return false;
     return _left < other._left ? ( other._left <= _right) : (_left <= other._right);
}

bool GenomicInterval::isContainedIn(const GenomicInterval &other, bool nonStrandness) const
{
     if( other._seq_id != _seq_id) return false;
     if( !nonStrandness && other._strand != Strand_t::StrandUnknown && _strand != Strand_t::StrandUnknown && other._strand != _strand) return false;
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
     if (!other.overlap(*this)) SError("this two interval does not overlap\n");
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
     if ( rhs._seq_id != _seq_id) return false;
     if( rhs._strand != Strand_t::StrandUnknown && _strand != Strand_t::StrandUnknown && rhs._strand != _strand) return false;
     return (_left == rhs._left && _right == rhs._right);
}

bool GenomicInterval::operator!=(const GenomicInterval& rhs) const
{
   return !(*this == rhs);
}
bool GenomicInterval::operator>(const GenomicInterval& rhs) const
{
     if ( rhs._seq_id != _seq_id) SError("cannot compare for different chrom\n");
     return (_left==rhs._left)?(_right>rhs._right):(_left>rhs._left);
}

bool GenomicInterval::operator<(const GenomicInterval& rhs) const
{
     if ( rhs._seq_id != _seq_id) SError("cannot compare for different chrom\n");
     return (_left == rhs._left)?(_right < rhs._right):(_left < rhs._left);
}

EmpDist::EmpDist(const std::vector<double>& pdf,
               const std::vector<double>& cdf,
               size_t mode_pos, double mean,
               double sd, size_t min, size_t max,
               FLD_source fld):
      _pdf(pdf), _cdf(cdf), _mode_pos(mode_pos),
      _mean(mean), _sd(sd), _min(min), _max(max),
      _fld_source(fld)
{}

void EmpDist::pdf(const std::vector<double>& pdf)
{
   _pdf = pdf;
}

double EmpDist::pdf(size_t i) const
{
   if( i >_max || i <_min)
      return 0.0;
   else
      return _pdf[i];
}

void EmpDist::cdf(const std::vector<double>& cdf)
{
   _cdf = cdf;
}

double EmpDist::cdf(size_t i) const
{
    if( i >_max || i <_min)
      return 0.0;
   else
      return _cdf[i];
}

void EmpDist::mode(size_t mode)
{
   _mode_pos = mode;
}

size_t EmpDist::mode() const
{
   return _mode_pos;
}


void EmpDist::mean(double mean){
   _mean = mean;
}

double EmpDist::mean() const
{
   return _mean;
}

void EmpDist::max(size_t max)
{
   _max = max;
}

size_t EmpDist::max() const
{
   return _max;
}

void EmpDist::min(size_t min)
{
   _min = min;
}
size_t EmpDist::min() const
{
   return _min;
}

void EmpDist::sd(double sd)
{
   _sd = sd;
}

double EmpDist::sd() const
{
   return _sd;
}

void EmpDist::fld_source(FLD_source fld){
   _fld_source = fld;
}

FLD_source EmpDist::fld_source() const
{
   return _fld_source;
}
