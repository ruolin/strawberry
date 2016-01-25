
/* **
 * Note:
 * Some functions are borrowed from Cufflinks: https://github.com/cole-trapnell-lab/cufflinks
 * **/

#include <libgen.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <errno.h>
#include "common.h"
bool WITH_BIAS_CORRECTION = false;
int kMaxGeneLength = 2500000;
int kMaxFragSpan = 1000000;
int kMaxFragPerCluster = 100000;
int kMaxIntronLength = 50000; // max-junction-splice-distance
int kMinIntronLength = 50; // min-junction-splice-distance
int kMaxIntronLen4ExtCluster = 3000; /*Do not extend the cluster if intron length large than this*/
int kMinTransLen = 200; //ignore isoforms if its length is too short.
int kMaxOlapDist = 30; // merge cluster if within this distance.
double kMaxSmallAnchor = 4;  // smallAnchor 4bp;
double kMinIsoformFrac = 0.05;
double kBinomialOverHangAlpha = 0.0;
float kMinJuncSupport = 1; // min number of spliced aligned reads for a valid intron
int kMinDist4ExonEdge = 5; // used in FlowNetwork::addWeight() for assigning
                                     // weights on non-intron edges.
double kMinDepth4Locus = 1; //used in ClusterFactory::finalizeCluster() to
                                        //select locus have enough reads covered.
double kMinDepth4Contig = 0;
int kMaxCoverGap1 = 200; // cover gap due the read depth.
int kMaxCoverGap2 = 50;
int kMaxReadNum4FD = 10000000;
//bool singleExon4FD = false;
//int kMinExonLen4FD = 200; // if singleExon4FD is used.
//int kMinExonCov4FD = 0;  // if singleExon4FD is used.
bool singleIso4FD = true;
double kInsertSizeMean = 0;
double bothStrandCutoff = 0.1; // ratio of reads from different strands. If larger than this value, it
                              // possibly indicates genes on different strands overlap.
double kInsertSizeSD = 0;
bool infer_the_other_end = true;
bool verbose = true;
bool kCombineShrotTransfrag = true;
std::string output_dir = ".";
std::string ref_gtf_filename = "";
std::string ref_fasta_file = "";
bool enforce_ref_models = false;
bool utilize_ref_models = false;
std::string tracking_log = "/tracking.log";
bool effective_len_norm = false;
float kIntronEdgeWeight = 2;
bool use_unique_hits = true;
bool use_paired_hits = false;
double standard_normal_cdf(double x)
/*
 * Implementation from http://www.johndcook.com/blog/cpp_phi/
 * */
{
   // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

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

const char* stripFileName(char *path)
{
   char *dummy = strdup(path);
   const char *dname = dirname(dummy);
   return dname;
};

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

void reverseString(char str[], int length)
{
    int start = 0;
    int end = length -1;
    while (start < end)
    {
        std::swap(*(str+start), *(str+end));
        start++;
        end--;
    }
}


/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */


int parseInt(const char* optarg, int lower, const char *errmsg, void (*print_help)()) {
    long l;
    char *endPtr= NULL;
    l = strtol(optarg, &endPtr, 10);
    if (endPtr != NULL) {
        if (l < lower) {
            std::cerr << errmsg << std::endl;
            print_help();
            exit(1);
        }
        return (int32_t)l;
    }
    std::cerr<< errmsg <<std::endl;
    print_help();
    exit(1);
    return -1;
}


float parseFloat(const char* optarg, float lower, float upper, const char *errmsg, void (*print_help)()) {
    float l;
    l = (float)atof(optarg);

    if (l < lower) {
        std::cerr << errmsg <<std::endl;
        print_help();
        exit(1);
    }

    if (l > upper)
    {
        std::cerr << errmsg << std::endl;
        print_help();
        exit(1);
    }

    return l;

    std::cerr << errmsg << std::endl;
    print_help();
    exit(1);
    return -1;
}

int mkpath(const char *s, mode_t mode)
{
   /*
    *  return code
    *  1: success
    *  0: already exist
    *  -1: cannot create
    */
    char *q, *r = NULL, *path = NULL, *up = NULL;
    int rv;

    rv = -1;
    if (strcmp(s, ".") == 0 || strcmp(s, "/") == 0)
        return (0);

    if ((path = strdup(s)) == NULL)
        exit(1);

    if ((q = strdup(s)) == NULL)
        exit(1);

    if ((r = dirname(q)) == NULL)
        goto out;

    if ((up = strdup(r)) == NULL)
        exit(1);

    if ((mkpath(up, mode) == -1) && (errno != EEXIST))
        goto out;

    if ((mkdir(path, mode) == -1) && (errno != EEXIST))
        rv = -1;
    else
        rv = 0;

out:
    if (up != NULL)
        free(up);
    free(q);
    free(path);
    return (rv);
}

// Implementation of itoa()
char* Sitoa(int num, char* str, int base)
{
    int i = 0;
    bool isNegative = false;

    /* Handle 0 explicitely, otherwise empty string is printed for 0 */
    if (num == 0)
    {
        str[i++] = '0';
        str[i] = '\0';
        return str;
    }

    // In standard itoa(), negative numbers are handled only with
    // base 10. Otherwise numbers are considered unsigned.
    if (num < 0 && base == 10)
    {
        isNegative = true;
        num = -num;
    }

    // Process individual digits
    while (num != 0)
    {
        int rem = num % base;
        str[i++] = (rem > 9)? (rem-10) + 'a' : rem + '0';
        num = num/base;
    }

    // If number is negative, append '-'
    if (isNegative)
        str[i++] = '-';

    str[i] = '\0'; // Append string terminator

    // Reverse the string
    reverseString(str, i);
    return str;
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
        delete buf;
        buf = new char[allocated];
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
     if (!other.overlap(*this)) {
        LOG_ERR("Calling overlapLen for two non-overlapping interval: ", _left, "-", _right,
              "\t", other._left, "-", other._right);
     }
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
     if ( rhs._seq_id != _seq_id) {
        return _seq_id > rhs._seq_id;
     }
     return (_left==rhs._left)?(_right>rhs._right):(_left>rhs._left);
}

bool GenomicInterval::operator<(const GenomicInterval& rhs) const
{
     if ( rhs._seq_id != _seq_id) {
        return _seq_id < rhs._seq_id;
     }
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

std::ostream& operator<<(std::ostream&os, const Strand_t& obj){
   switch(obj){
   case Strand_t::StrandPlus:
      os<<"+";
      break;
   case Strand_t::StrandMinus:
      os<<"-";
      break;
   default:
      os<<".";
   };
}
