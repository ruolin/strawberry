
 /*
 *  Created on: Jan 23, 2015
 *      Author: Ruolin Liu
 */

#ifndef FASTA_H_
#define FASTA_H_

#include<unordered_map>
#include<memory>

class FaRecord{
   /*
    * *  FASTA format: https://en.wikipedia.org/wiki/FASTA_format
 *
 *  >gi|31563518|ref|NP_852610.1| microtubule-associated proteins 1A/1B light chain 3A isoform b [Homo sapiens]
    MKMRFFSSPCGKAAVDPADRCKEVQQIRDQHPSKIPVIIERYKGEKQLPVLDKTKFLVPDHVNMSELVKI
    IRRRLQLNPTQAFFLLVNQHSMVSVSTPIADIYEQEKDEDGFLYMVYASQETFGFIRENE.
 *
 *  The sequence can be multiple lines.
 *
 *  One FaRecord is such an entry.
    */
public:
   std::string _seq_name; // e.g. gi|31563518|ref|NP_852610.1|
   uint _seq_len = 0;
   off_t _fpos = 0; // seq fpos in the file
   int _line_len = 0; //effective line length (without EoL)
   int _line_blen = 0; //length of line including EoL characters
   FaRecord(std::string name, uint len, off_t fpos, int llen, int blen):
      _seq_name(name),
      _seq_len(len),
      _fpos(fpos),
      _line_len(llen),
      _line_blen(blen){}
   FaRecord() = default;
};

class FaIndex{
   /*
    * Class represents a single fasta file.
    * In some cases, it is one fasta file per chromosome, therefore multiple objects.
    * In other cases, it is one fasta file contains all chromosomes, therefore one object.
    */
   std::string _fa_name; // fasta file name
   std::string _fai_name; // fasta index file name
   bool _haveFai;
public:
   std::unordered_map<std::string, FaRecord> _records; // map seq name to record.
   using FaRecord_p = std::unordered_map<std::string, FaRecord>::const_iterator;
   FaIndex(const char* fname, const char* finame=NULL);
   bool add_record(std::string seqname, const uint seqlen, const off_t fpos, const int linelen, const int lineblen);
   bool getRecord(const std::string& seqname, FaRecord &got) const;
   const std::string get_faidx_name() const;
   bool hasIndex();
   int loadIndex(); //return the number of record loaded
   int buildIndex(); //this function has not been implemented //return the number of record
   int writeIndex(); // return number of record which is stored
   int num_records() const;
};

class SubSeq{
public:
   uint _subseq_start;
   uint _subseq_len;
   char* _sequence;
   SubSeq():
      _subseq_start(1),
      _subseq_len(0),
      _sequence(nullptr){}
   ~SubSeq(){
      delete []_sequence;
      _sequence = nullptr;
   }
   SubSeq& operator=(const SubSeq &rhs) = delete;
   SubSeq(const SubSeq &other) = delete;
   SubSeq(SubSeq &&rhs):
      _subseq_start(rhs._subseq_start),
      _subseq_len(rhs._subseq_len),
      _sequence(rhs._sequence)
   {
      rhs._sequence = nullptr;
   }

   SubSeq& operator=(SubSeq &&rhs){
      if(this != &rhs){
         delete []_sequence;
         _subseq_start = rhs._subseq_start;
         _subseq_len = rhs._subseq_len;
         _sequence = rhs._sequence;
         rhs._sequence = nullptr;
      }
      return *this;
   }
   void setup(uint s, uint l);
};

class FaSeqGetter{
   /*
    * One seq at a time in memory
    */
   std::string _fname;
   FILE* _fh = nullptr;
   SubSeq _my_subseq;
   FaRecord _my_record;
public:
   const static int MAX_LEN_TO_READ = 0x20000000;
   FaSeqGetter() = default;
   void initiate(const std::string fname, const FaRecord &rec);
   std::string get_fname() const;
   // Default parameters mean loading a whole sequence
   // for the first time.
   uint loadSeq(uint start = 1, uint len = 0);
   std::string fetchSeq(uint start, uint len);
   ~FaSeqGetter();
   FaSeqGetter(const FaSeqGetter &other) = delete;
   FaSeqGetter& operator=(const FaSeqGetter &other) = delete;
   FaSeqGetter(FaSeqGetter &&other) = delete;
   FaSeqGetter& operator=(FaSeqGetter &&other) = delete;
};

class FaInterface{
   bool _has_load = false;
public:
   std::string _fa_path;
   std::unordered_map<std::string, std::unique_ptr<FaIndex>> _fa_indexes; //map fasta file name to faidx
   std::unordered_map<std::string, std::string> _seqname_2_fafile;   // map seq name to fasta file name.
   using ItFaidx = std::unordered_map<std::string, std::unique_ptr<FaIndex>>::iterator;
   void initiate(const char* fpath=nullptr);
   void load2FaSeqGetter(FaSeqGetter &getter, const std::string seqname);
   bool hasLoad() const{
      return _has_load;
   }
};


#endif /* FASTA_H_ */
