#include "genome.h"
GenomicInterval::GenomicInterval(){
      _start = 0;
      _end = 0;
      _chrom = NULL;
      _strand = 0;
   }
GenomicInterval::GenomicInterval(const char* chr, uint s, uint e, char o) {
      _chrom = Gstrdup(chr);
      _chrom = strlower(_chrom);
      if (s>e) { _start = e; _end = s;}
      else { _start = s; _end = e;}
      switch(o){
      case '+':
         _strand = strand_plus;
         break;
      case '-':
         _strand = strand_minus;
         break;
      case '.':
         _strand = strand_unknown;
         break;
      default:
         GError("Error: stand must be '+', '-' or '.'\n");
         break;
      }
   }

GenomicInterval::GenomicInterval( const GenomicInterval& other){
      _start = other._start;
      _end = other._end;
      _chrom = Gstrdup(other._chrom);
      _strand = other._strand;
   }


uint GenomicInterval::getStart(){
     return _start;
  }

uint GenomicInterval::getEnd(){
     return _end;
}

char GenomicInterval::getStrand(){
     return _strand;
}

char* GenomicInterval::getChrom(){
     return _chrom;
}

  //check for overlap with other segment
uint GenomicInterval::len() { return _end-_start+1; }



bool GenomicInterval::overlap(GenomicInterval* d, bool nonStrandness) {
     if( d == NULL) return false;
     if( strcmp(d->_chrom, _chrom)) return false;
     if( !nonStrandness && d->_strand != strand_unknown && _strand != strand_unknown && d->_strand != _strand) return false;
     return _start < d->_start ? ( d->_start <= _end) : (_start <= d->_end);
}

bool GenomicInterval::isContainedIn(GenomicInterval *d, bool nonStrandness){
     if( d == NULL) return false;
     if( strcmp(d->_chrom, _chrom)) return false;
     if( !nonStrandness && d->_strand != strand_unknown && _strand != strand_unknown && d->_strand != _strand) return false;
     if (_start < d->_start || _end > d->_end) return false;
     return true;
}

bool GenomicInterval::contain(GenomicInterval *d, bool nonStrandness){
     if (d == NULL) return false;
     return d->isContainedIn(this, nonStrandness);
}

  //return the length of overlap between two segments
uint GenomicInterval::overlapLen(GenomicInterval* r) {
     if (!r->overlap(this)) GError("this two interval does not overlap");
     if (_start<r->_start) {
        if (r->_start>_end) return 0;
        return (r->_end>_end) ? _end-r->_start+1 : r->_end-r->_start+1;
        }
       else { //r->start<=start
        if (_start>r->_end) return 0;
        return (r->_end<_end)? r->_end-_start+1 : _end-_start+1;
        }
}

GffAttr::GffAttr(int an_id, const char* av ): _attr_id(an_id){
   _attr_val = NULL;
   setAttrValue(av);
}
GffAttr:: ~GffAttr(){
   GFREE(_attr_val);
}

void GffAttr::setAttrValue(const char*av){
   if (_attr_val!=NULL) {
        GFREE(_attr_val);
        }
     if (av==NULL || av[0]==0) return;
     //trim spaces
     const char* vstart=av;
     while (*vstart==' ') vstart++;

     const char* vend=vstart;
     bool keep_dq=false;
     while (vend[1]!=0) {
        if (*vend==' ' && vend[1]!=' ') keep_dq=true;
        else if (*vend==';') keep_dq=true;
        vend++;
     }
     //remove spaces at the end:
     while (*vend==' ' && vend!=vstart) vend--;
     //practical clean-up: if it doesn't have any internal spaces just strip those useless double quotes
     if (!keep_dq && *vstart=='"' && *vend=='"') {
               vend--;
               vstart++;
     }
     _attr_val=Gstrdup(vstart, vend);
}

GffExon::GffExon( GenomicInterval iv, int exnum, string tid,
         double score, char phase):
  _exon_iv(iv),
  _exon_num(exnum),
  _parent_transcript_id(tid),
  _score(score),
  _phase(phase)
   {} //constructor

GffExon::~GffExon(){};

GffTranscript::GffTranscript(GenomicInterval iv, string tid):
  _trans_iv(iv), _trans_id(tid){
}
bool GffTranscript::setParentGeneID(string gid){
   if (_parent_gene_id.empty()){
   _parent_gene_id = gid;
   return true;
   }
   else{
      return false;
   }
}

void GffTranscript::addExon(exonPtr e){
   _exons.push_back(e);
}


GffGene::GffGene(GenomicInterval iv, string gid): _gene_iv(iv), _gene_id(gid){
}

void GffGene::addTranscript(transPtr &t){
	_transcripts.insert({t->_trans_id, t});
   _transcript_num++;
}

void GffLine::discardParent(){
   _num_parents=0;
   _parents_len=0;
  _parents=NULL;
}

char* GffLine::extractAttr(const char* attr, bool caseStrict, bool enforce_GTF2) {
  // This function is modified based on cufflinks2.0.0 gff.cpp
  //parse a key attribute and remove it from the info string
 //(only works for attributes that have values following them after ' ' or '=')
 static const char GTF2_ERR[]="Error parsing attribute %s ('\"' required) at GTF line:\n%s\n";
 int attrlen=strlen(attr);
 char cend=attr[attrlen-1];
 //char* pos = (caseStrict) ? strstr(info, attr) : strifind(info, attr);
 //must make sure attr is not found in quoted text
 char* pos=_info;
 char prevch=0;
 bool in_str=false;
 bool notfound=true;
 //int (*strcmpfn)(const char*, const char*, int) = caseStrict ? Gstrcmp : Gstricmp;
 while (notfound && *pos) {
   char ch=*pos;
   if (ch=='"') {
     in_str=!in_str;
     pos++;
     prevch=ch;
     continue;
     }
   char* temp = Gstrdup(pos);
   temp[attrlen]=0;
   if (!in_str && (prevch==0 || prevch==' ' || prevch == ';')
          && Gstricmp(attr, temp)==0) {
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
      }
   //not a match or in_str
   prevch=ch;
   pos++;
   }
 if (notfound) return NULL;
 char* vp=pos+attrlen;
 while (*vp==' ') vp++;
 if (*vp==';' || *vp==0)
      GError("Error parsing value of GFF attribute \"%s\", line:\n%s\n", attr, _dupline);
 bool dq_enclosed=false; //value string enclosed by double quotes
 if (*vp=='"') {
     dq_enclosed=true;
     vp++;
     }
 if (enforce_GTF2 && !dq_enclosed)
      GError(GTF2_ERR,attr, _dupline);
 char* vend=vp;
 if (dq_enclosed) {
    while (*vend!='"' && *vend!=';' && *vend!=0) vend++;
    }
 else {
    while (*vend!=';' && *vend!=0) vend++;
    }
 if (enforce_GTF2 && *vend!='"')
     GError(GTF2_ERR, attr, _dupline);
 char *r=Gstrdup(vp, vend-1);
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
   GMALLOC(_line,_llen+1);
   memcpy(_line, l,_llen+1);
   GMALLOC(_dupline, _llen+1);
   memcpy(_dupline,l,_llen+1);
   _skip=false;
   _chrom=NULL;
   _track=NULL;
   _ftype=NULL;
   _info=NULL;
   _parents=NULL;
   _num_parents=0;
   _parents_len=0;
   _is_gff3=false;
   //_is_cds=false; //for future
   _is_transcript=false;
   _is_exon=false;
   _is_gene=false;
   _exontype=0;
   _ID=NULL;
   _name=NULL;
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
   _chrom=t[0];
   _track=t[1];
   _ftype=t[2];
   _info=t[8];
   char* p=t[3];
   if(!parseUInt(p,_start)){
      GMessage("Warning: invalid start coordinate at line:\n%s\n",l);
      return;
   }
   p=t[4];
   if (!parseUInt(p,_end)){
      GMessage("Warning: invalid end coordinate at line:\n%s\n",l);
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
      if(!parseDouble(p,_score))
         GMessage("Warning: invalid feature score at line:\n%s\n",l);
         return;
   }
   _strand = *t[6];
   if(_strand!=strand_plus && _strand!=strand_minus
         && _strand!=strand_unknown){
      GMessage("Warning: parsing strand (%c) from GFF line:\n%s\n",_strand,l);
   }
   _phase=*t[7];
   _ID=NULL;
   strncpy(fnamelc, _ftype, 127);
   fnamelc[127]=0;
   strlower(fnamelc);
   bool is_t_data=false;
   if(strstr(fnamelc,"utr")!=NULL){
      _exontype=exGffUTR;
      _is_exon = false;
      is_t_data = true;
   }
   else if(strstr(fnamelc,"exon")!=NULL){
      _exontype = exGffExon;
      _is_exon = true;
      is_t_data = true;
   }
   else if (strstr(fnamelc, "stop") &&
         (strstr(fnamelc, "codon") || strstr(fnamelc, "cds"))){
      _exontype=exGffStop;
      _is_exon = false;
      is_t_data=true;
   }
   else if (strstr(fnamelc, "start") &&
         ((strstr(fnamelc, "codon")!=NULL) || strstr(fnamelc, "cds")!=NULL)){
      _exontype=exGffStart;
      _is_exon=false;
      is_t_data=true;
   }
   else if (strcmp(fnamelc, "cds")==0) {
      _exontype=exGffCDS;
      _is_exon=false;
      is_t_data=true;
   }
   else if (strstr(fnamelc, "gene")!=NULL) {
      _is_gene=true;
      is_t_data=true; //because its name will be attached to parented transcripts
   }
   else if (strstr(fnamelc,"rna")!=NULL || strstr(fnamelc,"transcript")!=NULL) {
      _is_transcript=true;
      is_t_data=true;
   }
   /*if(!is_t_data){
      return;
   }*/
   _ID=extractAttr("ID=",true);

   char* Parent=extractAttr("Parent=",true);
   _is_gff3=(_ID!=NULL || Parent!=NULL);
   if (_is_gff3) {
      if (_ID!=NULL) {
         //has ID attr so it's likely to be a parent feature
         //look for explicit gene name
         _name=extractAttr("name=");

         /* deal with messy name convention in gff3. We dont need for now
         if(_name==NULL){
            _name=extractAttr("gene_name=");
            if (_name==NULL) {
                _name=extractAttr("geneName=");
                if (_name==NULL) {
                    _name=extractAttr("gene_sym=");
                    if (_name==NULL) {
                            _name=extractAttr("gene=");
                    }
                }
            }
         }
         */

         if(_is_gene) GFREE(Parent); //TMI, we really don't care about gene Parents?
      } // has ID field
      if(Parent!=NULL){
         _num_parents=1;
         p=Parent;
         int last_delim_pos=-1;
         while(*p != ';' && *p!=0){
            if(*p ==',' && *(p+1)!=0 && *(p+1)!=';'){
               _num_parents++;
               last_delim_pos=(p-Parent);
            }
            p++;
         }
         _parents_len=p-Parent+1;
         GMALLOC(_parents,_num_parents*sizeof(char*));
         _parents[0]=Parent;
         int i=1;
         if(last_delim_pos>0){
            for (p=Parent+1;p<=Parent+last_delim_pos;++p){
               if(*p==','){ // substitute comma with string ending '\0'
                  char* ep=p-1;
                  while(*ep==' ' && ep>Parent) ep--;
                  *(ep+1)=0;
                  _parents[i]=p+1;
               }
            }
         }
      } // has Parent field
   } //GFF3
}

/*
GffLine::GffLine(GffLine *l){
   memcpy((void*)this, (void*)l, sizeof(GffLine));
   _line=NULL;
   GMALLOC(_line, _llen+1);
   memcpy(_line, l->_line, _llen+1);
   GMALLOC(_dupline, _llen+1);
   memcpy(_dupline, l->_dupline, _llen+1);
   _chrom=_line+(l->_chrom-l->_line);
   _track=_line+(l->_track-l->_line);
   _ftype=_line+(l->_ftype-l->_line);
   _info=_line+(l->_info-l->_line);
   if (l->_parents_len>0) {
      _parents_len=l->_parents_len;
      GMALLOC(_parents, _parents_len);
      memcpy(_parents, l->_parents, _parents_len);
   }
   if (l->_name!=NULL)
      _name=Gstrdup(l->_name);
   if (l->_ID!=NULL)
      _ID=Gstrdup(l->_ID);
}
*/

GffLine::GffLine(){
   _line=NULL;
   _dupline=NULL;
   _chrom=NULL;
   _track=NULL;
   _ftype=NULL;
   _start=0;
   _end=0;
   _info=NULL;
   _parents=NULL;
   _parents_len=0;
   _ID = NULL;
   _name = NULL;
   _strand = strand_unknown;
   _skip = false;
   _exontype=0;
   //_is_cds=false; // for future
   _is_gff3=false;
   _is_transcript=false;
   _is_gene=false;
   _is_exon=false;
   _llen=0;
   _phase=0;
   _score=0.0;
   _num_parents=0;
}

GffLine::~GffLine() {
     GFREE(_dupline);
     GFREE(_line);
     GFREE(_parents);
     GFREE(_ID);
     GFREE(_name);
     //GFREE(_chrom);
     //GFREE(_track);
     //GFREE(_ftype);
     //GFREE(_info);
    }

GffReader::GffReader(char* fname){
   _fname=Gstrdup(fname);
   _fh = fopen(fname,"rb");
   GMALLOC(_linebuf, GFF_LINELEN);
   _buflen=GFF_LINELEN-1;
   _fpos=0;
}

GffReader::~GffReader(){
   _fpos=0;
   GFREE(_linebuf);
   GFREE(_fname);
   fclose(_fh);
}

bool GffReader::nextGffLine(){
   if (_gffline) {
      //debug error
      GError("readAll() should free gffline after processing");
   }
   while(!_gffline){
      int llen=0;
      _buflen=GFF_LINELEN-1;
      char *l = fgetline(_linebuf, _buflen, _fh, &_fpos, &llen);
      if (l == NULL) return false; //end of file
      int ns=0; // first nonspace position
      while (l[ns]!=0 && isspace(l[ns])) ns++;
      if(l[ns]=='#' || llen<10) continue;
      _gffline = make_shared<GffLine>(l);
      if(_gffline->_skip){
         _gffline.reset();
         continue;
      }
      if(_gffline->_ID==NULL && _gffline->_parents==NULL){
         GMessage("Warning: malformed GFF line, no parent or record Id skipping\n");
         _gffline.reset();
      }

   }
   if(_gffline) return true;
   else return false;
}

void GffReader::readAll(){
   int exon_count=0;
   while(nextGffLine()){
      if(_gffline->_is_gene){
         exon_count = 0;
         //if(last_gene){
        // }
         GenomicInterval iv(_gffline->_chrom, _gffline->_start, _gffline->_end, _gffline->_strand);
         geneSharedPtr current_gene (new GffGene (iv, _gffline->_ID));
         _genes.insert({current_gene->_gene_id, current_gene});
      }
      else if(_gffline->_is_transcript){
         GenomicInterval iv(_gffline->_chrom, _gffline->_start, _gffline->_end, _gffline->_strand);
         transPtr current_trans (new GffTranscript(iv,_gffline->_ID));
         _genes[_gffline->_parents[0]]->addTranscript(current_trans);
         _trans2gene.insert({_gffline->_ID, _gffline->_parents[0]});
      }
      else if (_gffline->_is_exon){
         ++exon_count;
         GenomicInterval iv(_gffline->_chrom, _gffline->_start, _gffline->_end, _gffline->_strand);
         exonPtr current_exon(new GffExon(iv, exon_count, _gffline->_parents[0], _gffline->_score, _gffline->_phase));
         string geneID = _trans2gene[_gffline->_parents[0]];
         _genes[geneID]->_transcripts[_gffline->_parents[0]]->addExon(current_exon);
      }
      _gffline.reset();
   }
}


