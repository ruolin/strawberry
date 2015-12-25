/*
 * common.h
 *
 *  Created on: Jan 25, 2015
 *      Author: ruolin
 */

#ifndef COMMON_H_
#define COMMON_H_
#include<string>
#include<string.h>
#include<stdarg.h>
#include<vector>
#include <sys/stat.h>
#include <memory>
#include "logger.hpp"
typedef void* pointer;
typedef uint64_t ReadID;
typedef int RefID;

const static int kMaxInnerDist = 5000;
//const static int kMinTransLen =
const static int kMaxIntronLength = 30000;
const static double kSmallOverHangProp = 3/100.0;
const static double kMinIsoformFrac = 0.05;
const static double kBinomialOverHangAlpha = 0.0;
const static bool enforce_ref_models = false;
const static float kMinJuncSupport = 0.1; // min number of spliced aligned reads for a valid intron
const static int kMinDist4ExonEdge = 1; // used in FlowNetwork::addWeight() for assigning
                                        // weights on non-intron edges.
const static double kMinDepth4Locus = 0.1; //used in ClusterFactory::finalizeCluster() to
                                           //select locus have enough reads covered.
const static int kMaxCoverGap = 500;
const static int kMaxReadNum4FD = 5000;
const static double kInsertSizeMean = 200;
const static double kInsertSizeSD = 80;
const static bool infer_the_other_end = true;
const static bool verbose = true;
//const static int kMinExonLen = 5;
#define SFREE(ptr)       SFree((pointer*)(&ptr))

double standard_normal_cdf(double x);

int fileExists(const char* fname);

inline int64_t fileSize(const char* fpath)
{
   struct stat results;
   if (stat(fpath, &results) == 0)
      // The size of the file in bytes
      return (int64_t)results.st_size;
   else
      return 0;
}

bool endsWith (std::string const &fullString, std::string const &ending);

bool endsWith(const char* s, const char* suffix);

void split(const std::string& s, const std::string& delims, std::vector<std::string>& result);

inline void SFree(pointer* ptr)
{
   if (*ptr) free(*ptr);
   *ptr=NULL;
}
void reverseString(char str[], int length);
char* Sitoa(int num, char* str, int base);

inline void str2lower(char * str) {//changes string in place
  if (str==NULL) return;
  int i=0;
  while (str[i]!=0) { str[i]=tolower(str[i]); i++; }
}

inline void str2lower(std::string &str){
   if (str.empty()) return ;
   for(size_t i = 0; i<str.length(); i++)
      str[i] = tolower(str[i]);
}

inline void str2upper(char * str) {//changes string in place
  if (str==NULL) return;
  int i=0;
  while (str[i]!=0) { str[i]=toupper(str[i]); i++; }
}

inline void str2lupper(std::string &str){
   if (str.empty()) return ;
   for(size_t i = 0; i<str.length(); i++)
      str[i] = toupper(str[i]);
}

inline double getMedian(const vector<float> &vec){
   vector<float> dup = vec;
   double median = 0.0;
   size_t n = vec.size();
   if(n % 2 ==0)
      median = (dup[n/2] + dup[n/2-1]) / 2.0;
   else
      median = dup[n/2];
   return median;
}

int stricmp(const char* a, const char* b, int n);


bool SRealloc(pointer* ptr,unsigned long size);

inline bool overlaps_locally(uint lhs_left, uint lhs_right, uint rhs_left, uint rhs_right){

   if (lhs_left >= rhs_left && lhs_left < rhs_right)
      return true;
   if (lhs_right > rhs_left && lhs_right < rhs_right)
      return true;
   if (rhs_left >= lhs_left && rhs_left < lhs_right)
      return true;
   if (rhs_right > lhs_left && rhs_right < lhs_right)
      return true;
   return false;

}


//--------------------------------------------------------
// ************** simple line reading class for text files
class SlineReader {
protected:
   bool closeFile;
   int len; // number of characters read in a line
   int allocated;
   char* buf;
   bool isEOF;
   FILE* file;
   off_t filepos; //current position
   bool pushed; //pushed back
   int lcount; //line counter (read lines)
public:
   virtual char* chars() { return buf; }
   virtual char* line() { return buf; }
   virtual int readcount() { return lcount; } //number of lines read
   virtual void setFile(FILE* stream) { file=stream; }
   virtual int length() { return len; }
   virtual int size() { return len; } //same as size();
   virtual bool isEof() {return isEOF; }
   virtual bool eof() { return isEOF; }
   virtual off_t getfpos() { return filepos; }
   virtual off_t getFpos() { return filepos; }
   virtual char* nextLine() { return getLine(); }
   virtual char* getLine() { if (pushed) { pushed=false; return buf; }
                            else return getLine(file);  }
   virtual char* getLine(FILE* stream) {
                 if (pushed) { pushed=false; return buf; }
                          else return getLine(stream, filepos); }
   virtual char* getLine(FILE* stream, off_t& f_pos); //read a line from a stream and update
                           // the given file position
   virtual void pushBack() { if (lcount>0) pushed=true; } // "undo" the last getLine request
            // so the next call will in fact return the same line
   SlineReader(const char* fname) {
      FILE* f=fopen(fname, "rb");
      if (f==NULL) {
         LOG_ERR("Error opening file: ",fname);
      }
      closeFile=true;
      init(f);
      }
   SlineReader(FILE* stream=NULL, off_t fpos=0) {
     closeFile=false;
     init(stream,fpos);
     }
   virtual void init(FILE* stream, off_t fpos=0) {
     len=0;
     isEOF=false;
     allocated=1024;
     buf = new char[allocated];
     lcount=0;
     buf[0]=0;
     file=stream;
     filepos=fpos;
     pushed=false;
     }
   virtual ~SlineReader() {
     delete[] buf;
     if (closeFile) fclose(file);
     }
};

enum class Strand_t: char{
   StrandUnknown,
   StrandPlus,
   StrandMinus,
   StrandBoth
};

std::ostream& operator<<(std::ostream&os, const Strand_t& obj);

class GenomicInterval {
private:
   uint _left = 1; // left < right always!
   uint _right = 0;
   int _seq_id = -1; // _chrom is seq id. -1 is for unmapped read
   Strand_t _strand;

public:
   GenomicInterval()=default;
   GenomicInterval(int chr,
                  uint l,
                  uint r,
                  Strand_t o);

  uint left() const;
  uint right() const;
  void set_seq_id(int id);
  void set_left(uint l);
  void set_right(uint r);
  Strand_t strand() const;
  int seq_id() const;
  uint len() const;

  //check for overlap with other segment
  bool overlap(const GenomicInterval &other, bool nonStrandness = true) const;

  bool isContainedIn(const GenomicInterval &other, bool nonStrandness = true) const;

  bool contain(const GenomicInterval &other, bool nonStrandness = true) const;

  //return the length of overlap between two segments
  uint overlapLen(const GenomicInterval& other) const;
  //comparison operators required for sorting
  bool operator==(const GenomicInterval& d) const;
  bool operator!=(const GenomicInterval& rhs) const;
  bool operator>(const GenomicInterval& d) const;
  bool operator<(const GenomicInterval& d) const;
};

enum FLD_source
{
    LEARNED,
    USER,
    DEFAULT
};

class EmpDist{
   // truncated dist between min and max
   std::vector<double> _pdf;
   std::vector<double> _cdf;
   size_t _mode_pos;
   double _mean;
   double _sd;
   size_t _min;
   size_t _max;
   FLD_source _fld_source;
public:
   EmpDist(const std::vector<double>& pdf, const std::vector<double>& cdf,
         size_t mode_pos, double mean, double sd, size_t min, size_t max, FLD_source fld);
   EmpDist() = default;

   void pdf(const std::vector<double>& pdf);
   double pdf(size_t i) const;

   void cdf(const std::vector<double>& cdf);
   double cdf(size_t i) const;

   void mode(size_t mode);
   size_t mode() const;

   void mean(double mean);
   double mean() const;

   void max(size_t max);
   size_t max() const;

   void min(size_t min);
   size_t min() const;

   void sd(double sd);
   double sd() const;

   void fld_source(FLD_source fld);
   FLD_source fld_source() const;
};

enum strandedness_t
{
   UNKNOWN_STRANDNESS,
   STRANDED,
   UNSTRANDED
};

enum mate_strand_orien_t
{
    UNKNOWN_MATE_ORIENTATION,
    MATES_POINT_TOWARD,
    MATES_POINT_SAME,
    MATES_POINT_AWAY,
    UNPAIRED,
};

enum mate_strand_mapping_t
{
   //should be either FR or RF or both
   FF,
   FR,
   RF,
   RR
};

enum platform_t
{
    UNKNOWN_PLATFORM,
    ILLUMINA,
    SOLID
};

struct AssayProperties
{
   platform_t _platform;
   strandedness_t _strandedness;
   mate_strand_mapping_t _mate_strand_mapping;
   mate_strand_orien_t _mate_strand_orien;
   long double _total_mapped_mass;
   long double _norm_mapped_mass;
   std::unique_ptr<EmpDist> _frag_len_dist;
   std::string _condition_name;
   std::string _file_path;
   int _num_replicates;
   AssayProperties() = default;
};

template <class ForwardIterator, class OutputIterator>
  ForwardIterator unique2 (ForwardIterator first, ForwardIterator last, OutputIterator out)
{
  if (first==last) return last;
  typedef typename std::iterator_traits<ForwardIterator>::value_type value;
  ForwardIterator result = first;
  ForwardIterator begin = first;
  ForwardIterator end =first;
  while (++first != last)
  {
    if (!(*result == *first)){
       while(begin != end){
            value val = *(++begin);
            *(out++) = val;
       }
      *(++result)=*first;
      begin = first;
      end = first;
    }
    else{
       end = first;
    }
  }
   while(begin != end){
            value val = *(++begin);
            *(out++) = val;
   }
  return ++result;
}



struct IntronTable{
   uint left;
   uint right;
   float total_junc_reads;
   float small_span_read;
   double median_depth;
   //vector<float> doc;
   IntronTable(uint l, uint r):
      left(l),
      right(r),
      total_junc_reads(0.0),
      small_span_read(0.0),
      median_depth(0)
   {}
   bool operator==(const IntronTable & rhs){
      return (left == rhs.left && right == rhs.right);
   }
   bool operator<(const IntronTable &rhs){
      if(left != rhs.left)
         return left < rhs.left;
      if(right != rhs.right)
         return right < rhs.right;
      return false;
   }
   static bool overlap(const IntronTable& lhs, const IntronTable& rhs){
      return overlaps_locally(lhs.left, lhs.right, rhs.left, rhs.right);
   }
};



#endif /* COMMON_H_ */
