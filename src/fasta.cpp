#include "fasta.h"
#include "common.h"
#include <dirent.h>
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <mutex>

using namespace std;
// for now this function free whole sequence and setup a
// new _my_subseq object for FaSeqGetter. This is not an efficient way.
void SubSeq::setup(uint s, uint l){
    if (_sequence){ // if not empty
       delete [] _sequence;
    };
   _subseq_start = s;
   _subseq_len = l;
   _sequence = new char[l+1];
}

//initialize .fa file and .fai file name
FaIndex::FaIndex(const char* fname, const char* finame){
   if(fileExists(fname) != 2) {
      LOG_ERR("Error: fasta file ", fname, " not found");
      exit(1);
   }
   if (fileSize(fname)<=0) {
      LOG_ERR("Error: invalid fasta file",fname);
      exit(1);
   }
   _fa_name.assign(fname);
   if(finame) {
      _fai_name.assign(finame);
      if(fileExists(finame) == 2 and fileSize(finame) > 0){
         loadIndex();
      } else{
         buildIndex();
      }
   }
}

int FaIndex::num_records() const { return _records.size();}
int FaIndex::loadIndex(){
   _records.clear();
   _haveFai = false;
   FILE* fi = fopen(_fai_name.c_str(), "rb");
   if(!fi){
      LOG_ERR("cannot open fasta index file for reading ",_fai_name.c_str());
      LOG_ERR("Currently creating fasta index is not supported. Please use samtools to build fasta index.");
      exit(1);
   }
   SlineReader fl(fi);
   char* line = nullptr;
   while( (line = fl.nextLine()) != nullptr ){
      if(*line == '#') continue;
      size_t idx = strcspn(line, " \t"); // first occurrence of space or tab delimiter
      if(idx == strlen(line)) {
         LOG_ERR("Error parsing fasta index line: ",line);
         exit(1);
      }
      char *p = (char*)(line+idx); //s now holds the genomic sequence name
      *(p++) = 0;
      uint len=0;
      int line_len=0, bline_len=0;
      #ifdef __WIN32__
         long offset = -1;
         sscanf(p, "%d%ld%d%d", &len, &offset, &line_len, &bline_len);
      #else
         long long offset = -1;
         sscanf(p, "%d%lld%d%d", &len, &offset, &line_len, &bline_len);
      #endif
      if(len==0 || line_len==0 || bline_len==0 || line_len > bline_len){
         LOG_ERR("Error parsing fasta index line: ",line);
         exit(1);
      }
      #ifdef DEBUG
         //printf("%s\t%d\t%lld\t%d\t%d\n", line, len ,offset, line_len, bline_len);
      #endif
      //str2lower(line);
      add_record(line, len ,offset, line_len, bline_len);
   }
   fclose(fi);
   return _records.size();

}
// NEED to be done. Currently faIdx is generated by samtools.
int FaIndex::buildIndex(){
   return 0;
}

const string FaIndex::get_faidx_name() const {return _fai_name;}

bool FaIndex::add_record(string seqname, const uint seqlen, const off_t fpos, const int linelen, const int lineblen){
   pair<FaRecord_p, bool> res;
   res = _records.insert( make_pair(seqname, FaRecord(seqname, seqlen, fpos, linelen, lineblen) ) );
   if(!res.second){
      LOG_WARN("duplicated seqname ", seqname, " in fasta file");
   }

   return res.second;
}

bool FaIndex::getRecord(const string &seqname, FaRecord &got) const{
   auto it = _records.find(seqname);
   if(it !=  _records.end()){
      got = it->second;
      return true;
   }
   return false;
}

FaSeqGetter::~FaSeqGetter(){
   if(_fh != nullptr) fclose(_fh);
   _fh = nullptr;
}

void FaSeqGetter::initiate(const string fname, const FaRecord &rec)
{
   _fname = fname;
   _my_record = rec;
   if(_fh != nullptr){
      fclose(_fh);
      _fh == nullptr;
   }
   _fh = fopen(fname.c_str(),"rb");
}

string FaSeqGetter::get_fname() const {return _fname;}

uint FaSeqGetter::loadSeq(uint start, uint len){
   // get faIdx info from record
   //start is 1-based genomic coordinate within current fasta sequence
   int line_len = _my_record._line_len;
   int line_blen = _my_record._line_blen;
   int line_endlen = line_blen - line_len;
   off_t orig_start = _my_record._fpos;
   uint seq_len = _my_record._seq_len;
   if(seq_len == 0){
      LOG_ERR("Empty or zero-length fasta record ", _my_record._seq_name);
      exit(1);
   }
   uint start_line_number = (start-1) / line_len;
   int start_char_in_line = (start-1) % line_len;
   off_t f_start = orig_start + (start_line_number*line_blen) + start_char_in_line;
   uint toread = len;
   if(toread == 0) toread = seq_len - start +1;
   if(toread > MAX_LEN_TO_READ) toread = MAX_LEN_TO_READ;

   _my_subseq.setup(start, toread);
   //assert(start >= _my_subseq._subseq_start);
   char* cur_char_p = _my_subseq._sequence;
   fseeko(_fh, f_start, SEEK_SET);
   int act_read_len = 0;
   uint already_read_len = 0;
   if(start_char_in_line > 0){ // read first line (partially if indeed)
      int should_read_len = line_len - start_char_in_line;
      if (should_read_len > toread) should_read_len = toread; // in case we need just a few chars
      act_read_len = fread((void*)cur_char_p, 1, should_read_len, _fh);
      if( act_read_len < should_read_len){
         LOG_ERR("reading ",_fname, " encountered a premature eof. Please check input.");
         exit(1);
      }
      toread -= act_read_len;
      cur_char_p += act_read_len;
      already_read_len += act_read_len;
      fseeko(_fh, line_endlen, SEEK_CUR);
   }
   while(toread >= line_len){
      act_read_len = fread((void*)cur_char_p, 1, line_len, _fh);
      if(act_read_len < line_len){
         cerr<<"reading "<<_fname<<" encountered a premature eof. Please check input."<<endl;
         fclose(_fh);
         exit(1);
      }
      toread -= act_read_len;
      cur_char_p += act_read_len;
      already_read_len += act_read_len;
      fseeko(_fh, line_endlen, SEEK_CUR);
   }
   if(toread>0){ // read last line
      act_read_len = fread((void*)cur_char_p, 1, toread, _fh);
      if(act_read_len < toread){
         cerr<<"reading "<<_fname<<" encountered a premature eof. Please check input."<<endl;;
         fclose(_fh);
         exit(1);
      }
      already_read_len += act_read_len;
   }

   return already_read_len;
}

string FaSeqGetter::fetchSeq(uint start, uint len){
   char str[len+1];
   strncpy(str, _my_subseq._sequence+start-1, len);
   str[len] = 0;
   return string(str);
}


void FaInterface::initiate(const char* fpath){
   // fpath should be check non-empty in the caller
   _fa_path.assign(fpath);
   pair<ItFaidx, bool> ret;
   switch(fileExists(fpath)){
   case 0:
   {
      cerr<<"File or directory "<<fpath<<" does not exist!"<<endl;
      exit(1);
   }
   case 2: // is a file. One file multiple chromosomes.
   {
      if(endsWith(_fa_path,".fai") ){ // if .fai file is given instead of .fa
         string fa_file_name = _fa_path.substr(0, _fa_path.length()-4);
         if(!fileExists(fa_file_name.c_str() )  ){
            cerr<<"Cannot find fasta file for index file "<<fpath<<endl;
            exit(1);
         }else{
            ret = _fa_indexes.insert(make_pair(fa_file_name, unique_ptr<FaIndex> (new FaIndex(fa_file_name.c_str(), fpath) ) ) );
            assert(ret.second);
         }
      } else if(endsWith(_fa_path, ".fa") || endsWith(_fa_path, ".fasta")){
         string fai_name = _fa_path+".fai";
         ret = _fa_indexes.insert(make_pair(string(fpath), unique_ptr<FaIndex>(new FaIndex(fpath, fai_name.c_str()) ) )  );
         assert(ret.second);

      } else {
         cerr<<"Cannot find .fasta or .fa file"<<endl;
         exit(1);
      }

      // this for loop initialize _seqname_2_fafile object
      unique_ptr<FaIndex> &faidx_ptr = ret.first->second;
      const string &fasta_file_name  = ret.first->first;
      for(auto record = faidx_ptr->_records.begin(); record != faidx_ptr->_records.end(); ++record){
         pair<unordered_map<string, string>::iterator, bool> ret_it;
         ret_it = _seqname_2_fafile.insert(make_pair(record->first,fasta_file_name));
         if(!ret_it.second){
            LOG_ERR("Please checking fasta file ", fasta_file_name,  "for possible duplicated sequence names" );
         }
      }//end for loop

      break;
   }
   case 1:  // is a directory. one file one chromosome.
   {
      DIR *dir;
      struct dirent *ent;
      dir = opendir(fpath);
      while((ent = readdir(dir)) != NULL){
         if(endsWith(ent->d_name, ".fa") || endsWith(ent->d_name, ".fasta")){
            char fai_name[200];
            fai_name[0] = 0;
            strcat(fai_name, fpath);
            if(!endsWith(fpath, "/")) strcat(fai_name,"/");
            strcat(fai_name, ent->d_name);
            string fa_file_name(fai_name);
            strcat(fai_name, ".fai");
            if(strlen(fai_name) > 199) {
               cerr<<"file name is too long "<<fai_name<<endl;
               exit(0);
            }
            // if index file exists
            if( fileExists(fai_name) ==2 ) {
               ret = _fa_indexes.insert(make_pair(fa_file_name, unique_ptr<FaIndex>(new FaIndex(fa_file_name.c_str(), fai_name) ) )  );
               assert(ret.second);
               // this for loop initialize _seqname_2_fafile object
               unique_ptr<FaIndex> &faidx_ptr = ret.first->second;
               const string &fasta_file_name  = ret.first->first;
               for(auto record = faidx_ptr->_records.begin(); record != faidx_ptr->_records.end(); ++record){
                  pair<unordered_map<string, string>::iterator, bool> ret_it;
                  ret_it = _seqname_2_fafile.insert(make_pair(record->first, fasta_file_name));
                  if(!ret_it.second){
                     LOG_ERR("Please checking fasta file ", fasta_file_name,  "for possible duplicated sequence names" );
                  }
               }//end for loop
            }
            else{
               cerr<<"Warning: fasta file "<<fa_file_name<<" lack index file!"<<endl;
               if(!system(NULL)){
                  cerr<<"processor is not available"<<endl;
                  exit(0);
               }
               const char* samtools = "samtools faidx ";
               char* command = (char*) malloc(1 + strlen(samtools) + fa_file_name.size());
               strcpy(command, samtools);
               strcat(command, fa_file_name.c_str());
               cerr<<"Now using samtools to build fasta index file"<<endl;
               system(command);

               if(fileExists(fai_name) !=2){
                  cerr<<"Unable to call samtools. Please check you have completely install samtools!"<<endl;
               }

               ret = _fa_indexes.insert(make_pair(fa_file_name, unique_ptr<FaIndex>(new FaIndex(fa_file_name.c_str(), fai_name) ) )  );
               assert(ret.second);
               // this for loop initialize _seqname_2_fafile object
               unique_ptr<FaIndex> &faidx_ptr = ret.first->second;
               const string &fasta_file_name  = ret.first->first;
               for(auto record = faidx_ptr->_records.begin(); record != faidx_ptr->_records.end(); ++record){
                  pair<unordered_map<string, string>::iterator, bool> ret_it;
                  ret_it = _seqname_2_fafile.insert(make_pair(record->first, fasta_file_name));
                  if(!ret_it.second){
                     LOG_ERR("Please checking fasta file ", fasta_file_name,  "for possible duplicated sequence names" );
                  }
               }//end for loop
            }
         }
      }
      closedir(dir);
      break;
   }
   default:
      cerr<<"Error: not a valid file or directory "<<fpath<<endl;;
      break;
   }
   cerr<<"Load "<<_seqname_2_fafile.size()<<" reference fasta"<<endl;
   _has_load = true;

}

void FaInterface::load2FaSeqGetter(FaSeqGetter &getter, const string seqname){
   auto it_fa_file_name = _seqname_2_fafile.find(seqname);
   if(it_fa_file_name == _seqname_2_fafile.end()){
      cerr<<"Reference sequence name "<<seqname<<" cannot be found in fasta file. Please check fasta file header line."<<endl;
      //cerr<<_seqname_2_fafile.begin()->first<<endl;
      exit(0);
   }
   string fa_file_name = it_fa_file_name->second;
   auto it_faidx = _fa_indexes.find(fa_file_name);
   assert(it_faidx != _fa_indexes.end());
   FaRecord rec;
   if(it_faidx->second->getRecord(seqname, rec))
      getter.initiate(fa_file_name, rec);
   else{
      cerr<<"Fetching seq name "<<seqname<< " failed!"<<endl;
      exit(0);
   }
}
