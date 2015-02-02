#include "gff.h"
#include<cstring>
#include<cassert>
#include<algorithm>
#ifdef DEBUG
   #include<iostream>
   #include<stdio.h>
#endif



void GffLine::extractAttr(const string attr, string &val) {
   //parse a key attribute and remove it from the info string
   //(only works for attributes that have values following them after ' ' or '=')
   //static const char GTF2_ERR[]="Error parsing attribute %s ('\"' required) at GTF line:\n%s\n";
   int attrlen=attr.length();
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
      string temp(pos,attrlen);
      str2lower(temp);
      if (!in_str && (prevch==0 || prevch==' ' || prevch == ';') && attr == temp){ //attr match found

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
   if (notfound) return ;
   char* vp=pos+attrlen;
   while (*vp==' ') vp++;
   if (*vp==';' || *vp==0)
      SError("Error parsing value of GFF attribute \"%s\", line:\n%s\n", attr.c_str(), _dupline);
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

   val = string(vp, vend-vp);

   //-- now remove this attribute from the info string
   while (*vend!=0 && (*vend=='"' || *vend==';' || *vend==' ')) vend++;
   if (*vend==0) vend--;
   for (char *src=vend, *dest=pos;;src++,dest++) {
      *dest=*src;
      if (*src==0) break;
   }
   return;
}


static char fnamelc[128];
GffLine::GffLine(const char* l)
{
//   printf("enter in gffline constructor\n");
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
   extractAttr("id=", _ID);
   extractAttr("parent=", _parent);
   _is_gff3=(!_ID.empty() || !_parent.empty());
   if (_is_gff3) { // is gff3
      if (!_ID.empty()) { // has ID field
         //has ID attr so it's likely to be a parent feature
         //look for explicit gene name
         extractAttr("name=", _name);
         //deal with messy name convention in gff3. We dont need for now
         if(_name.empty()){
            extractAttr("gene_name=", _name);
            if (_name.empty()) {
                extractAttr("genename=", _name);
                if (_name.empty()) {
                    extractAttr("gene_sym=", _name);
                    if (_name.empty()) {
                            extractAttr("gene=", _name);
                    }
                }
            }
         }
      } // end has ID field
      if(!_parent.empty()){
         split(_parent, ",", _parents);
      }
   } // end is gff3

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

int GffInfoVec::addInfo(const string name){
   auto it = _name2id.find(name);
   if(it != _name2id.end()) return it->second;
   else{
      int idx = _gff_info_vec.size();
      GffInfo ginfo(idx, name);
      _gff_info_vec.push_back(ginfo);
      _name2id.insert(pair<string,int>(name, idx));
      return idx;
   }
   return -1;
}

unique_ptr<GffInfoTable> GffObj::_infotable = unique_ptr<GffInfoTable> (new GffInfoTable());

GffObj::GffObj(LinePtr gl, GffReader & greader):
   _iv(GenomicInterval(gl->_chrom, gl->_start, gl->_end, gl->_strand)),
   _score(gl->_score),
   _phase(gl->_phase),
   _source(gl->_source),
   _greader(greader)
{
   _seq_id = _infotable->_seq_names->addInfo(gl->_chrom);

}

GffLoci::GffLoci(LinePtr gl, GffReader & greader):
      GffObj(gl, greader),
      _gene_id (gl->_ID),
      _gene_name (gl->_name)
{}

GffmRNA::GffmRNA(LinePtr gl, GffReader & greader):
      GffObj(gl, greader),
      _transcript_id(gl->_ID),
      _transcript_name(gl->_name)
{
    if(gl->_parents.size() != 1)
       SMessage("No parent or multiple parents for a mRNA object in: %s\n", gl->_dupline);
    else
       _parent_gene = gl->_parents[0];
}

GffExon::GffExon(LinePtr gl, GffReader & greader):
      GffObj(gl, greader),
      _exon_id(gl->_ID),
      _exon_name(gl->_name)
{
   for(auto str : gl->_parents){
      _parent_mrnas.push_back(str);
   }
}

GffReader::GffReader(vector<unique_ptr<GffSeqData>> & gseqs, const char* fname):
      SlineReader(fname),
      _fname(string(fname)),
      _g_seqs(gseqs)
{}

bool GffReader::nextGffLine(){
   while(true){
      const char *l = nextLine();
      if (l == NULL) {

         return false; //end of file
      }
      int ns=0; // first nonspace position
      while (l[ns]!=0 && isspace(l[ns])) ns++;
      if(l[ns]=='#' ||len<10) continue;
      _gfline. reset (new GffLine(l));
      if(_gfline->_skip){
         continue;
      }
      if(_gfline->_ID.empty() && _gfline->_parent.empty()){
         SMessage("Warning: malformed GFF line, %s\n",_gfline->_dupline);
         continue;
      }
      break;
   }
   return true;
}

void GffReader::readAll(){
   string last_seq_name;
   int cur_gseq = 0;
   while(nextGffLine()){
      if( _gfline->_chrom != last_seq_name){
         last_seq_name = _gfline->_chrom;
         unique_ptr<GffSeqData> g_seq(new GffSeqData());
         _g_seqs.push_back(move(g_seq));
         cur_gseq = _g_seqs.size()-1;
      }
      switch(_gfline->_feat_type)
      {
      case GENE:
       {
         GffLoci loci(_gfline, *this);
         _g_seqs[cur_gseq]->_genes.push_back(move(loci));

         break;
       }
      case mRNA:
       {

         GffmRNA mrna(_gfline, *this);
         GffmRNA *cur_mrna = NULL;
         if(mrna.strand() == kStrandPlus){
            _g_seqs[cur_gseq]->_forward_rnas.push_back(mrna);
            cur_mrna = &_g_seqs[cur_gseq]->_forward_rnas.back();
         }
         else if(mrna.strand() == kStrandMinus){
            _g_seqs[cur_gseq] -> _reverse_rnas.push_back(mrna);
            cur_mrna = &_g_seqs[cur_gseq]->_reverse_rnas.back();
         }
         else{
            _g_seqs[cur_gseq] -> _unstranded_rnas.push_back(mrna);
            cur_mrna = &_g_seqs[cur_gseq]->_unstranded_rnas.back();
         }

         // it is most possible that the last gene is the parent.
         string p_gene = mrna.parent_gene();
         if( _g_seqs[cur_gseq]->last_gene()._gene_id == p_gene){
            _g_seqs[cur_gseq]->last_gene()._mrnas.push_back(cur_mrna);
         }
         else{ // In most cases this will not be executed. search for the gene list in gseq
            SMessage("gff file not in order. line: \n", _gfline->_dupline);
            auto it = find_if(_g_seqs[cur_gseq]->_genes.begin(), _g_seqs[cur_gseq]->_genes.end(),
                  [&p_gene](GffLoci const& g){return g._gene_id == p_gene;} );
            if(it == _g_seqs[cur_gseq]->_genes.end()) {
               SError("orphan mRNA found in GFF file: %s \n", _gfline->_dupline);
            } else {
               it->_mrnas.push_back(cur_mrna);
            }
         }
         break;
       }
      case EXON:
       {
          GffExon exon(_gfline, *this);
          GffExon *cur_exon = NULL;
       }
      default:
         break;
      }
   }
}


