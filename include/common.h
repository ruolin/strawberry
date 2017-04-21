/*
 * common.h
 *
 *  Created on: Jan 25, 2015
 *      Author: Ruolin Liu
 */

#ifndef COMMON_H_
#define COMMON_H_
#include<string>
#include<string.h>
#include<stdarg.h>
#include<vector>
#include <sys/stat.h>
#include <memory>
#include <cmath>
#include <iomanip>
#include "logger.hpp"
typedef void* pointer;
typedef uint64_t ReadID;
typedef int RefID;

extern bool SINGLE_END_EXP;
extern bool BIAS_CORRECTION;
extern int kMaxGeneLength;
extern int kMaxFragSpan;
extern int kMaxFragPerCluster;
extern int kMaxIntronLength; // max-junction-splice-distance
extern int kMinIntronLength; // min-junction-splice-distance
extern int kMaxIntronLen4ExtCluster; /*Do not extend the cluster if intron length large than this*/
extern unsigned int kMinExonLen;
extern int kMinTransLen; //ignore isoforms if its length is too short.
extern int kMaxOlapDist; // merge cluster if within this distance.
extern double kMaxSmallAnchor;  // smallAnchor 4bp;
extern double kMinIsoformFrac;
extern double kBinomialOverHangAlpha;
extern bool enforce_ref_models;
extern bool utilize_ref_models;
extern bool no_assembly;
extern float kMinJuncSupport; // min number of spliced aligned reads for a valid intron
extern int kMinDist4ExonEdge; // used in FlowNetwork::addWeight() for assigning
                                        // weights on non-intron edges.
extern double kMinDepth4Locus; //used in ClusterFactory::finalizeCluster() to
                                           //select locus have enough reads covered.
extern double kMinDepth4Contig;
extern double kMinExonDoc;
extern int kMaxCoverGap1;
extern int kMaxCoverGap2;
extern int kMaxReadNum4FD;
//extern int kMinExonLen4FD;
//extern int kMinExonCov4FD;
//extern bool singleExon4FD;

extern int num_threads;

extern double kInsertSizeMean;
extern double kInsertSizeSD;
extern bool infer_the_other_end;
extern bool verbose;
extern double kSimDepthCorrect;
extern bool kCombineShrotTransfrag;
extern double bothStrandCutoff;
extern std::string output_dir;
extern std::string ref_gtf_filename;
extern std::string ref_fasta_file;
extern std::string tracking_log;
extern std::string frag_context_out;
extern bool print_frag_context;
extern bool effective_len_norm;
extern float kIntronEdgeWeight;
extern bool use_only_unique_hits;
//extern bool use_only_paired_hits;
extern bool use_threads;
extern bool filter_by_expression;
extern bool weight_bias;
#define SFREE(ptr)       SFree((pointer*)(&ptr))
#define ENABLE_THREADS 1

double standard_normal_cdf(double x);

template <typename T>
T normal_pdf(T x, T m, T s)
{
    static const T inv_sqrt_2pi = 0.3989422804014327;
    T a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-T(0.5) * a * a);
}


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

double gc_content(const std::string& str);
void split(const std::string& s, const std::string& delims, std::vector<std::string>& result);

inline std::string fileName(const std::string& s) {
   std::vector<std::string> fields;
   split(s, "/", fields);
   return fields.back();
}

const char* stripFileName(char *path);

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

inline double getMedian(const std::vector<float> &vec){
   std::vector<float> dup = vec;
   double median = 0.0;
   size_t n = vec.size();
   if(n % 2 ==0)
      median = (dup[n/2] + dup[n/2-1]) / 2.0;
   else
      median = dup[n/2];
   return median;
}

int stricmp(const char* a, const char* b, int n);

inline bool overlaps_locally(uint lhs_left, uint lhs_right, uint rhs_left, uint rhs_right)
{
   return lhs_left <= rhs_right && rhs_left <= lhs_right;
}

int parseInt(const char* optarg,
          int lower,
          const char *errmsg,
          void (*print_help)());

float parseFloat(const char* optarg,
             float lower,
             float upper,
             const char *errmsg,
             void (*print_help)());
int mkpath(const char *s, mode_t mode);

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

template<typename str>
std::ostream& operator<<(std::ostream& os, const std::vector<str>& vec) {
   os<<"[";
   bool flag = true;
   for (const auto& item: vec) {
      if (flag) {
         os << item;
         flag = false;
      }
      else {
         os<<", "<<item;
      }
   }
   os<<"]"<<std::endl;
   return os;
}

inline void pretty_print(FILE* file, const std::vector<std::string>& vec, const std::string sep) {
   std::string out_str;
   bool flag = true;
   for (const auto& item: vec) {
      if (flag) {
         out_str += item;
         flag = false;
      }
      else {
         out_str += sep;
         out_str += item;
      }
   }
   out_str += "\n";
   fprintf(file, out_str.c_str());
}

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
   std::ostringstream out;
   out << std::setprecision(n) << a_value;
   return out.str();
}

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
