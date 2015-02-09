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
      if (!in_str && (prevch==0 || prevch==' ' || prevch == ';')
            && stricmp(attr.c_str(),pos, attrlen)==0){ //attr match found

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
   _line[_llen] = 0;
   memcpy(_line, l,_llen+1);
   _dupline =  new char[_llen+1];
   _dupline[_llen] = 0;
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
   if(_strand!=GenomicInterval::kStrandPlus && _strand!=GenomicInterval::kStrandMinus
         && _strand!=GenomicInterval::kStrandUnknown){
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
   else if (strstr(fnamelc,"rna")!=NULL || strstr(fnamelc,"transcript")!=NULL) {
      _feat_type = mRNA;
   }
   else if (strstr(fnamelc, "gene")!=NULL) {
      _feat_type = GENE;
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
   _strand = GenomicInterval::kStrandUnknown;
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
   // the id start from 0;
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
   _score(gl->_score),
   _phase(gl->_phase),
   _source(gl->_source),
   _greader(greader)
{
   int _seq_id = _infotable->_seq_names->addInfo(gl->_chrom);
   _iv = GenomicInterval(_seq_id, gl->_start, gl->_end, gl->_strand);
}

GffLoci::GffLoci(LinePtr gl, GffReader & greader):
      GffObj(gl, greader),
      _gene_id (gl->_ID),
      _gene_name (gl->_name)
{}

GffmRNA* GffLoci::getRNA(const string rna) {
   // in mast cases, we only need to look at the last mrna in the vector.
   if(_mrnas.back()->_transcript_id == rna){
      return _mrnas.back();
   } else{
      for(auto it = _mrnas.begin(); it != _mrnas.end(); ++it){
         if((*it)->_transcript_id == rna) return *it;
      }
      SMessage("Gene that contains %s does not show up in front of %s in gff file\n", rna.c_str(), rna.c_str());
   }
   return NULL;
}

void GffLoci::add_exon(exonPtr exon){
   if(num_mRNAs() == 1){
      if(_non_dup_exons.empty()){
         _mrnas[0]->add_exon( &(*exon));
         _non_dup_exons.push_back(move(exon));
      }
      else{
#ifdef DEBUG
            assert( *_non_dup_exons.back() < *exon);
#endif
      _non_dup_exons.push_back(move(exon));
      _mrnas[0]->add_exon( &(*exon));
      }
   }
   else{
      auto it = lower_bound(_non_dup_exons.begin(), _non_dup_exons.end(), exon,\
            [&](const exonPtr &lhs, const exonPtr &rhs){return *lhs < *rhs;});
      if( it == _non_dup_exons.cend()){
         _non_dup_exons.insert(it, move(exon));
      }else{
         if((**it) == *exon){
            (*it)->add_parent_mrna(exon->_parent_mrnas);
         }
         else{
            _non_dup_exons.insert(it, move(exon));
         }
      }
   }
}

GffmRNA::GffmRNA(LinePtr gl, GffLoci* gene, GffReader & greader):
      GffObj(gl, greader),
      _parent(gene),
      _transcript_id(gl->_ID),
      _transcript_name(gl->_name)
{}

GffExon::GffExon(LinePtr gl, GffmRNA* mrna, GffLoci* const gene, GffReader & greader):
      GffObj(gl, greader),
      _parent_gene(gene),
      _exon_id(gl->_ID),
      _exon_name(gl->_name)

{
   _parent_mrnas.push_back(mrna);
}

GffLoci* GffSeqData::findGene(const string gene_id){
   if( _genes.back()->_gene_id == gene_id){
      return &(*_genes.back());
   } else{
      for(auto it = _genes.begin(); it!= _genes.end(); it++){
         if( (*it)->_gene_id == gene_id) return &(*(*it));
      }
      SError("Gff file does not contain gene_id: %s \n", gene_id.c_str());
      return NULL;
   }
}

GffmRNA* GffSeqData::findmRNA(const string mrna_id, const char strand){
   switch(strand)
   {
      case '+':
      {
         if(_forward_rnas.back()->_transcript_id == mrna_id){
            return &(*_forward_rnas.back());
         } else {
            for(auto it = _forward_rnas.begin(); it!= _forward_rnas.end(); it++){
               if( (*it)->_transcript_id == mrna_id) return &(*(*it));
            }
            SError("GFF error: exon's parent %s cannot be found!\n", mrna_id.c_str());
         }
         break;
      }
      case '-':
      {
         if(_reverse_rnas.back()->_transcript_id == mrna_id){
            return &(*_reverse_rnas.back());
         } else {
            for(auto it = _reverse_rnas.begin(); it!= _reverse_rnas.end(); it++){
               if( (*it)->_transcript_id == mrna_id) return &(*(*it));
            }
            SError("GFF error: exon's parent %s cannot be found!\n", mrna_id.c_str());
         }
         break;
      }
      default:
      {
         if(_unstranded_rnas.back()->_transcript_id == mrna_id){
            return &(*_unstranded_rnas.back());
         } else {
            for(auto it = _unstranded_rnas.begin(); it!= _unstranded_rnas.end(); it++){
               if( (*it)->_transcript_id == mrna_id) return &(*(*it));
            }
            SError("GFF error: exon's parent %s cannot be found!\n", mrna_id.c_str());
         }
         break;
      }
   }
   return NULL;
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
   GffSeqData* gseq = nullptr;
   while(nextGffLine()){
      if( _gfline->_chrom != last_seq_name){
         last_seq_name = _gfline->_chrom;
         unique_ptr<GffSeqData> g_seq(new GffSeqData(last_seq_name));
         GffReader::addGseq(move(g_seq));
         gseq = &(*_g_seqs.back());
      }
      switch(_gfline->_feat_type)
      {
      case GENE:
       {
          unique_ptr<GffLoci> gene (new GffLoci(_gfline, *this));
          //SError("first seq id: %d \n", gene->get_seq_id());
         gseq->addGene(move(gene));
         break;
       }
      case mRNA:
       {

         GffLoci* gene = gseq->findGene(_gfline->parent());

         unique_ptr<GffmRNA> mrna( new GffmRNA(_gfline, gene, *this));
         GffmRNA *cur_mrna = NULL;
         if(mrna->strand() == GenomicInterval::kStrandPlus){
            gseq->addPlusRNA(move(mrna));
            cur_mrna = gseq->last_f_rna();
         } else if(mrna->strand() == GenomicInterval::kStrandMinus){
            gseq->addMinusRNA(move(mrna));
            cur_mrna = gseq->last_r_rna();
         } else{
            gseq->addUnstrandedRNA(move(mrna));
            cur_mrna = gseq->last_u_rna();
         }
         // mrna now is released
         // it is most possible that the last gene is the parent.
         gene->add_mRNA(cur_mrna);
         break;
       }
      case EXON:
       {
          GffmRNA* mrna = gseq->findmRNA(_gfline->parent(), _gfline->_strand);
          GffLoci *const gene  = mrna->getParentGene();
          exonPtr exon (new GffExon(_gfline, mrna, gene , *this));

          gene->add_exon(move(exon));
          //GffExon exon(_gfline, _g_seqs[cur_gseq]->last_gene(), *this);
          //GffExon *cur_exon = NULL;
       }
      default:
         break;
      }
   }
}


