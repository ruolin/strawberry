#include "gff.h"
#include<cstring>
#include "common.h"




char* GffLine::extractAttr(const char* attr) {
  // This function is modified based on cufflinks2.0.0 gff.cpp
  //parse a key attribute and remove it from the info string
 //(only works for attributes that have values following them after ' ' or '=')
 static const char GTF2_ERR[]="Error parsing attribute %s ('\"' required) at GTF line:\n%s\n";
 int attrlen=strlen(attr);
 char cend=attr[attrlen-1];
 //must make sure attr is not found in quoted text
 char* pos=_info;
 char prevch=0;
 bool in_str=false;
 bool notfound=true;
 while (notfound && *pos) {
   char ch=*pos;
   if (ch=='"') {
     in_str=!in_str;
     pos++;
     prevch=ch;
     continue;
     }
   char* temp = nullptr;
   strcpy(temp, pos);
   temp[attrlen]=0;
   str2lower(temp);
   if (!in_str && (prevch==0 || prevch==' ' || prevch == ';')
          && strcmp(attr, temp)==0) {
      //attr match found
      //check for word boundary on right
      char* epos=pos+attrlen;
      if (cend=='=' || cend==' ' || *epos==0 || *epos==' ') {
        notfound=false;
        break;
        }
      //not a perfect match, move on
      pos=epos;
      prevch=*(pos-1);
      continue;
      }//not a match or in_str
   prevch=ch;
   pos++;
   }
 if (notfound) return NULL;
 char* vp=pos+attrlen;
 while (*vp==' ') vp++;
 if (*vp==';' || *vp==0)
      SError("Error parsing value of GFF attribute \"%s\", line:\n%s\n", attr, _dupline);
 bool dq_enclosed=false; //value string enclosed by double quotes
 if (*vp=='"') {
     dq_enclosed=true;
     vp++;
     }

 char* vend=vp;
 if (dq_enclosed) {
    while (*vend!='"' && *vend!=';' && *vend!=0) vend++;
    }
 else {
    while (*vend!=';' && *vend!=0) vend++;
    }

 char *r = NULL;
 strncpy(r, vp, vend-vp);
 //-- now remove this attribute from the info string
 while (*vend!=0 && (*vend=='"' || *vend==';' || *vend==' ')) vend++;
 if (*vend==0) vend--;
 for (char *src=vend, *dest=pos;;src++,dest++) {
   *dest=*src;
   if (*src==0) break;
   }
 return r;
}


static char fnamelc[128];
GffLine::GffLine(const char* l) {
   _llen=strlen(l);
   _line = new char[_llen+1];
   memcpy(_line, l,_llen+1);
   _dupline =  new char[_llen+1];
   memcpy(_dupline,l,_llen+1);
   _skip=false;
   _is_gff3=false;
   _info = nullptr;
   _feat_type = OTHERS;
   //_is_cds=false; //for future
   _start=0;
   _end=0;
   char *t[9];
   int i=0;
   int tidx=1;
   t[0] = _line;
   while(_line[i]!=0){
      if(_line[i]=='\t'){
         _line[i]=0;
         t[tidx]=_line+i+1;
         tidx++;
         if(tidx>8) break;
      }
      i++;
   }
   if(tidx<8){
      return;
   }
   _chrom=string(t[0]);
   _source=string(t[1]);
   _gffline_type = string(t[2]);
   _info=t[8];
   char* p=t[3];
   _start = (uint) atol(p);
   if(_start == 0){
      SMessage("Warning: invalid start coordinate at line:\n%s\n",l);
      return;
   }
   p=t[4];
   _end = (uint) atol(p);
   if (_end == 0){
      SMessage("Warning: invalid end coordinate at line:\n%s\n",l);
      return;
   }
   if (_end<_start) {
      uint tem = _end;
      _end = _start;
      _start = tem;
   }; //make sure start>=end, always
   p=t[5];
   if (p[0]=='.' && p[1]==0) {
    _score=0;
    }
   else{
      _score = atof(p);
      if(_score == 0.0)
         SMessage("Warning: invalid feature score at line:\n%s\n",l);
         return;
   }
   _strand = *t[6];
   if(_strand!=kStrandPlus && _strand!=kStrandMinus
         && _strand!=kStrandUnknown){
      SMessage("Warning: parsing strand (%c) from GFF line:\n%s\n",_strand,l);
   }
   _phase=*t[7];
   strncpy(fnamelc, t[2], 127);
   fnamelc[127]=0;
   str2lower(fnamelc);
   if(strstr(fnamelc,"utr")!=NULL){
      _feat_type = UTR;
   }
   else if(strstr(fnamelc,"exon")!=NULL){
      _feat_type = EXON;
   }
   else if (strstr(fnamelc, "stop") &&
         (strstr(fnamelc, "codon") || strstr(fnamelc, "cds"))){
      _feat_type = STOP_CODON;
   }
   else if (strstr(fnamelc, "start") &&
         ((strstr(fnamelc, "codon")!=NULL) || strstr(fnamelc, "cds")!=NULL)){
      _feat_type = START_CODON;
   }
   else if (strcmp(fnamelc, "cds")==0) {
      _feat_type = CDS;
   }
   else if (strstr(fnamelc, "gene")!=NULL) {
      _feat_type = GENE;
   }
   else if (strstr(fnamelc,"rna")!=NULL || strstr(fnamelc,"transcript")!=NULL) {
      _feat_type = mRNA;
   }

   char *id = extractAttr("id=");
   char *parent = extractAttr("parent=");
   char *name = nullptr;
   _is_gff3=(id!=NULL ||parent!=NULL);
   if (_is_gff3) {
      if (id!=NULL) {
         //has ID attr so it's likely to be a parent feature
         //look for explicit gene name
         name=extractAttr("name=");
         //deal with messy name convention in gff3. We dont need for now
         if(name==NULL){
            name=extractAttr("gene_name=");
            if (name==NULL) {
                name=extractAttr("genename=");
                if (name==NULL) {
                    name=extractAttr("gene_sym=");
                    if (name==NULL) {
                            name=extractAttr("gene=");
                    }
                }
            }
         }
         _name = string(name);
         free(name);
         name = nullptr;
         _ID = string(id);
         free(id);
         id = nullptr;
         if(_feat_type == GENE){
            free(parent); //we really don't care about gene Parents?
            parent = NULL;
         }
      } // has ID field
      if(parent!=NULL){
         _parent = string(parent);
         free(parent);
         parent = NULL;
      } // has Parent field
   } //GFF3
}



GffLine::GffLine(){
   _line=NULL;
   _dupline=NULL;
   _start=0;
   _end=0;
   _info=NULL;
   _strand = kStrandUnknown;
   _skip = false;
   //_is_cds=false; // for future
   _is_gff3=false;
   _llen=0;
   _phase=0;
   _score=0.0;
}

GffLine::~GffLine() {
   delete[] _line;
   delete[] _dupline;
   _line = NULL;
   _dupline = NULL;
   _info = NULL;
   }

GffLine::GffLine(const GffLine &other):
      _llen (other._llen),
      _end (other._end),
      _start (other._start),
      _strand (other._strand),
      _skip (other._skip),
      _score (other._score),
      _chrom(other._chrom),
      _gffline_type(other._gffline_type),
      _feat_type (other._feat_type),
      _phase (other._phase),
      _ID (other._ID),
      _name(other._name),
      _parent(other._parent),
      _parents(other._parents)
{
   _line = new char[_llen+1];
   strncpy(_line, other._line, _llen);
   _line[_llen] = 0;
   _dupline = new char[_llen+1];
   strncpy(_dupline, other._line, _llen);
   _dupline[_llen] = 0;
   _info = other._info;
}

GffLine::GffLine(GffLine &&other):
      _llen (other._llen),
      _end (other._end),
      _start (other._start),
      _strand (other._strand),
      _skip (other._skip),
      _score (other._score),
      _chrom(move(other._chrom)),
      _gffline_type(move(other._gffline_type)),
      _feat_type (other._feat_type),
      _phase (other._phase),
      _ID (move(other._ID)),
      _name(move(other._name)),
      _parent(move(other._parent)),
      _parents(move(other._parents))
{
   _line = other._line;
   _dupline = other._dupline;
   _info = other._info;
   other._line = NULL;
   other._dupline = NULL;
   other._info = NULL;
}

GffReader::GffReader(const char* fname):
      SlineReader(fname),
      _fname(string(fname))
{}

GffReader::~GffReader(){
   _fpos=0;
   fclose(_fh);
}

bool GffReader::nextGffLine(){
   assert (!_gffline);
   while(!_gffline){
      int llen = 0;
      char *l = nextLine();
      if (l == NULL) return false; //end of file
      int ns=0; // first nonspace position
      while (l[ns]!=0 && isspace(l[ns])) ns++;
      if(l[ns]=='#' ||llen<10) continue;
      _gffline = make_shared<GffLine>(l);
      if(_gffline->_skip){
         _gffline.reset();
         continue;
      }
      if(_gffline->_ID==NULL && _gffline->_parents==NULL){
         Message("Warning: malformed GFF line, no parent or record Id skipping\n");
         _gffline.reset();
         continue;
      }
   }
   if(_gffline) return true;
   else return false;
}

void GffReader::readAll(){
   sl(_fh, _fpos);
   int exon_count=0;
   while(nextGffLine()){
      _gffline.reset();
   }
}


