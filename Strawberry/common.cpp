/*
 * common.cpp
 *
 *  Created on: Feb 2, 2015
 *      Author: ruolin
 */
#include "common.h"

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
