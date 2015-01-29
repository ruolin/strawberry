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
#include <sys/stat.h>

static char msg[4069];
typedef void* pointer;
#define SMALLOC(ptr,size)  if (!SMalloc((pointer*)(&ptr),size)) \
                                     SError("Error allocating memory.\n")
#define SFREE(ptr)       SFree((pointer*)(&ptr))
#define SREALLOC(ptr,size) if (!SRealloc((pointer*)(&ptr),size)) \
                                     SError("Error allocating memory.\n")


inline int fileExists(const char* fname) {
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

inline int64_t fileSize(const char* fpath) {
  struct stat results;
  if (stat(fpath, &results) == 0)
      // The size of the file in bytes
      return (int64_t)results.st_size;
  else
    return 0;
}

inline bool endsWith (string const &fullString, string const &ending)
{
   if (fullString.length() >= ending.length()) {
      return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
   } else {
      return false;
   }
}

inline bool endsWith(const char* s, const char* suffix) {
 if (suffix==NULL || s==NULL) return false;
 if (suffix[0]==0) return true; //special case: empty suffix
 int j=strlen(suffix)-1;
 int i=strlen(s)-1;
 if (i<j) return false;
 while (j>=0 && s[i]==suffix[j]) { i--; j--; }
 return (j==-1);
 }


void SFree(pointer* ptr){
  if (*ptr) free(*ptr);
  *ptr=NULL;
  }

inline void SError(const char* format,...){
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

inline void SMessage(const char* format,...){
  va_list arguments;
  va_start(arguments,format);
  vsprintf(msg,format,arguments);
  va_end(arguments);
  #ifdef __WIN32__
    OutputDebugString(msg);
  #endif
  fprintf(stderr,"%s",msg);fflush(stderr);
  }

inline bool SMalloc(pointer* ptr,unsigned long size){
  //GASSERT(ptr);
  if (size!=0) *ptr=malloc(size);
  return *ptr!=NULL;
  }


inline bool SRealloc(pointer* ptr,unsigned long size){
  //GASSERT(ptr);
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
   }//malloc
  else {//realloc
   void *p=realloc(*ptr,size);
   if (p) {
       *ptr=p;
       return true;
       }
   return false;
   }
 }


//--------------------------------------------------------
// ************** simple line reading class for text files
class SlineReader {
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
   char* chars() { return buf; }
   char* line() { return buf; }
   int readcount() { return lcount; } //number of lines read
   void setFile(FILE* stream) { file=stream; }
   int length() { return len; }
   int size() { return len; } //same as size();
   bool isEof() {return isEOF; }
   bool eof() { return isEOF; }
   off_t getfpos() { return filepos; }
   off_t getFpos() { return filepos; }
   char* nextLine() { return getLine(); }
   char* getLine() { if (pushed) { pushed=false; return buf; }
                            else return getLine(file);  }
   char* getLine(FILE* stream) {
                 if (pushed) { pushed=false; return buf; }
                          else return getLine(stream, filepos); }
   char* getLine(FILE* stream, off_t& f_pos); //read a line from a stream and update
                           // the given file position
   void pushBack() { if (lcount>0) pushed=true; } // "undo" the last getLine request
            // so the next call will in fact return the same line
   SlineReader(const char* fname) {
      FILE* f=fopen(fname, "rb");
      if (f==NULL) SError("Error opening file '%s'!\n",fname);
      closeFile=true;
      init(f);
      }
   SlineReader(FILE* stream=NULL, off_t fpos=0) {
     closeFile=false;
     init(stream,fpos);
     }
   void init(FILE* stream, off_t fpos=0) {
     len=0;
     isEOF=false;
     allocated=1024;
     SMALLOC(buf,allocated);
     lcount=0;
     buf[0]=0;
     file=stream;
     filepos=fpos;
     pushed=false;
     }
   ~SlineReader() {
     SFREE(buf);
     if (closeFile) fclose(file);
     }
};

char* SlineReader::getLine(FILE* stream, off_t& f_pos) {
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


#endif /* COMMON_H_ */
