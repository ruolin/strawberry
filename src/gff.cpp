#include "gff.h"
#include<cstring>
#include<cassert>
#include<algorithm>
#ifdef DEBUG
   #include<iostream>
   #include<stdio.h>
#endif

using namespace std;
unique_ptr<GffInfoTable> GffObj::_infotable = unique_ptr<GffInfoTable> (new GffInfoTable());

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
      std::cerr<<"Error parsing value of GFF attribute "<<attr<< "on line "<< _dupline << std::endl;
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


static char feat_type[128];
GffLine::GffLine(const char* l)
{
//   printf("enter in gffline constructor\n");
   //printf("%s\n", l);
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
   str2lower(_chrom);
   _source=string(t[1]);
   _gffline_type = string(t[2]);
   _info=t[8];
   char* p=t[3];
   _start = (uint) atol(p);
   if(_start == 0){
      std::cerr<<"invalid start coordinate at line:\n"<<l << std::endl;
      return;
   }
   p=t[4];
   _end = (uint) atol(p);
   if (_end == 0){
      std::cerr <<"invalid end coordinate at line:\n"<<l << std::endl;
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
         std::cerr<<"invalid feature score at line:\n"<<l << std::endl;
         return;
   }
   switch(*t[6]){
      case '+':
         _strand = Strand_t::StrandPlus; break;
      case '-':
         _strand = Strand_t::StrandMinus; break;
      default:
         _strand = Strand_t::StrandUnknown; break;
   }

   _phase=*t[7];
   strncpy(feat_type, t[2], 127);
   feat_type[127]=0;
   str2lower(feat_type);
   if(strstr(feat_type,"utr")!=NULL){
      _feat_type = UTR;
   }
   else if(strstr(feat_type,"exon")!=NULL){
      _feat_type = EXON;
   }
   else if (strstr(feat_type, "stop") &&
         (strstr(feat_type, "codon") || strstr(feat_type, "cds"))){
      _feat_type = STOP_CODON;
   }
   else if (strstr(feat_type, "start") &&
         ((strstr(feat_type, "codon")!=NULL) || strstr(feat_type, "cds")!=NULL)){
      _feat_type = START_CODON;
   }
   else if (strcmp(feat_type, "cds")==0) {
      _feat_type = CDS;
   }
   else if (strstr(feat_type,"rna")!=NULL ||
            strstr(feat_type,"transcript")!=NULL ||
            strstr(feat_type, "RNA") != NULL) {
      _feat_type = mRNA;
   }
   else if (strstr(feat_type, "gene")!=NULL) {
      _feat_type = GENE;
   }
   else {
      _feat_type = OTHERS;
      _skip = true;
   }
   extractAttr("id=", _ID);
   extractAttr("parent=", _parent);
   _is_gff3=(!_ID.empty() || !_parent.empty());
   if (_is_gff3) { // is gff3
      if (!_ID.empty()) { // has ID field
         if (_ID.find("transcript") != std::string::npos) _feat_type = mRNA;
         else if (_ID.find("gene") != std::string::npos) _feat_type = GENE;
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
   else{
      extractAttr("gene_id", _gene_id);
      extractAttr("gene_name", _gene_name);
      extractAttr("transcript_id", _transcript_id);
      // crazy hack
      if (_feat_type == GENE) {
         _ID = _gene_id;
         _name = _gene_name;
      }
      else {
         _ID = _transcript_id;
      }
   }
}



GffLine::GffLine(){
   _line=NULL;
   _dupline=NULL;
   _start=0;
   _end=0;
   _info=NULL;
   _strand = Strand_t::StrandUnknown;
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

GffObj::GffObj(LinePtr gl /*,GffReader & greader*/):
   _score(gl->_score),
   _phase(gl->_phase),
   _source(gl->_source)
   //_greader(greader)
{
   uint _seq_id = _infotable->_seq_names->addInfo(gl->_chrom);
   _iv = GenomicInterval(_seq_id, gl->_start, gl->_end, gl->_strand);
}


GffmRNA* GffLoci::getRNA(const string rna) {
   // in mast cases, we only need to look at the last mrna in the vector.
   if(_mrnas.back()->_transcript_id == rna){
      return _mrnas.back();
   } else{
      for(auto it = _mrnas.begin(); it != _mrnas.end(); ++it){
         if((*it)->_transcript_id == rna) return *it;
      }
      std::cerr<<"Can not find the parent of mRNA "<< rna << std::endl;
   }
   return NULL;
}

void GffLoci::add_exon(exonPtr exon, GffmRNA* exon_parent_mrna){
   if(num_mRNAs() == 1){
      if(_non_dup_exons.empty()){
         exon_parent_mrna->add_exon( &(*exon));
         assert((*exon)._iv.left());
         _non_dup_exons.push_back(move(exon));
      }
      else{
#ifdef DEBUG
            assert( *_non_dup_exons.back() < *exon);
#endif
      exon_parent_mrna->add_exon( &(*exon));
      _non_dup_exons.push_back(move(exon));
      }
   }
   else{
      auto it = lower_bound(_non_dup_exons.begin(), _non_dup_exons.end(), exon,\
            [&](const exonPtr &lhs, const exonPtr &rhs){return *lhs < *rhs;});
      if( it == _non_dup_exons.cend()){
         exon_parent_mrna->add_exon( &(*exon) );
         _non_dup_exons.insert(it, move(exon));
      } else{
         if((**it) == *exon){
            exon_parent_mrna->add_exon( &(**it));
            (*it)->add_parent_mrna(exon->_parent_mrnas);
         }
         else{
            exon_parent_mrna->add_exon(&(*exon));
            _non_dup_exons.insert(it, move(exon));
         }
      }
   }
}


GffExon::GffExon(LinePtr gl, GffmRNA* mrna, GffLoci* const gene /*, GffReader & greader*/):
      GffObj(gl /*, greader*/),
      _parent_gene(gene),
      //_exon_id(gl->_ID),
      _exon_name(gl->_name)

{
   _parent_mrnas.push_back(mrna);
}

GffLoci* GffTree::findGene(const string gene_id){
   if (_genes.empty()) return nullptr;
   if( _genes.back()->_gene_id == gene_id){
      return &(*_genes.back());
   } else{
      for(auto it = _genes.begin(); it!= _genes.end(); it++){
         if( (*it)->_gene_id == gene_id) return &(*(*it));
      }
      //unique_ptr<GffLoci> gene (new GffLoci(gene_id));
      //addGene(move(gene));
      return nullptr;
   }
}

GffmRNA* GffTree::findmRNA(const string mrna_id, const Strand_t strand){
   switch(strand)
   {
      case Strand_t::StrandPlus:
      {
         if(_forward_rnas.empty()) return nullptr;
         if(_forward_rnas.back()->_transcript_id == mrna_id){
            return &(*_forward_rnas.back());
         } else {
            for(auto it = _forward_rnas.begin(); it!= _forward_rnas.end(); it++){
               if( (*it)->_transcript_id == mrna_id) return &(*(*it));
            }
         }
         break;
      }
      case Strand_t::StrandMinus:
      {
         if(_reverse_rnas.empty()) return nullptr;
         if(_reverse_rnas.back()->_transcript_id == mrna_id){
            return &(*_reverse_rnas.back());
         } else {
            for(auto it = _reverse_rnas.begin(); it!= _reverse_rnas.end(); it++){
               if( (*it)->_transcript_id == mrna_id) return &(*(*it));
            }
         }
         break;
      }
      case Strand_t::StrandUnknown:
      {
         if(_unstranded_rnas.empty()) return nullptr;
         if(_unstranded_rnas.back()->_transcript_id == mrna_id){
            return &(*_unstranded_rnas.back());
         } else {
            for(auto it = _unstranded_rnas.begin(); it!= _unstranded_rnas.end(); it++){
               if( (*it)->_transcript_id == mrna_id) return &(*(*it));
            }
         }
         break;
      }
      //default:
         //assert(false);
   }
   return nullptr;
}

GffReader::GffReader(const char* fname, FILE* stream):
      SlineReader(stream), _fname(string(fname))
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
      _gfline.reset(new GffLine(l));
      if(_gfline->_skip){
         continue;
      }
      break;
   }
   return true;
}

void GffReader::readAll(){
   string previous_chrom;
   GffTree* gseq = nullptr;
   while(nextGffLine()){
      if( _gfline->_chrom != previous_chrom){
         previous_chrom = _gfline->_chrom;
         unique_ptr<GffTree> g_seq(new GffTree(previous_chrom));
         GffReader::addGseq(move(g_seq));
         gseq = &(*_g_seqs.back());
      }
      switch(_gfline->_feat_type)
      {
      case GENE:
       {
         unique_ptr<GffLoci> gene (new GffLoci(_gfline));
          //SError("first seq id: %d \n", gene->get_seq_id());
         gseq->addGene(move(gene));
         break;
       }
      case mRNA:
       {

         GffLoci* gene = gseq->findGene(_gfline->parent());
         if (gene == nullptr) break;

         unique_ptr<GffmRNA> mrna( new GffmRNA(_gfline, gene));
         GffmRNA *cur_mrna = NULL;
         if(mrna->strand() == Strand_t::StrandPlus){
            gseq->addPlusRNA(move(mrna));
            cur_mrna = gseq->last_f_rna();
         } else if(mrna->strand() == Strand_t::StrandMinus){
            gseq->addMinusRNA(move(mrna));
            cur_mrna = gseq->last_r_rna();
         } else if(mrna->strand() == Strand_t::StrandUnknown){
            gseq->addUnstrandedRNA(move(mrna));
            cur_mrna = gseq->last_u_rna();
         }
         else{
            assert(false);
         }
         // mrna now is released
         // it is most possible that the last gene is the parent.
         gene->add_mRNA(cur_mrna);
         break;
       }
      case EXON:
       {
          string parent = _gfline->parent();
          if (parent.empty()) { // gtf format
             parent = _gfline->_transcript_id;
             //cerr<<parent<<" parent\n";
             //cerr<<_gfline->_gene_id<<" gene id\n";
          }
          assert(!parent.empty());
          GffmRNA* mrna = gseq->findmRNA(parent, _gfline->_strand);
          GffLoci* gene = nullptr;
          if (mrna == nullptr) {
             gene = gseq->findGene(_gfline->_gene_id);
             if (gene == nullptr) {
               unique_ptr<GffLoci> locus(new GffLoci(_gfline->_gene_id));
               gseq->addGene(move(locus));
               gene = gseq->findGene(_gfline->_gene_id);
               assert(gene != nullptr);
             }
             mrnaPtr transcript(new GffmRNA(_gfline, gene));
             if(transcript->strand() == Strand_t::StrandPlus){
                gseq->addPlusRNA(move(transcript));
                mrna = gseq->last_f_rna();
             } else if(transcript->strand() == Strand_t::StrandMinus){
                gseq->addMinusRNA(move(transcript));
                mrna = gseq->last_r_rna();
             } else if(transcript->strand() == Strand_t::StrandUnknown){
                gseq->addUnstrandedRNA(move(transcript));
                mrna = gseq->last_u_rna();
             }

          } else {
             gene = mrna->getParentGene();
          }
          //cerr<<mrna->_transcript_id<<" transcript id"<<endl;
          exonPtr exon(new GffExon(_gfline, mrna, gene));
          assert(exon->_iv.left());
          gene->add_exon(move(exon), mrna);
          //GffExon exon(_gfline, _g_seqs[cur_gseq]->last_gene(), *this);
          //GffExon *cur_exon = NULL;
       }
      default:
         break;
      }
   }
//#ifdef DEBUG
//   cout<< gseq->_genes[0]->_non_dup_exons.size()<<endl;
//   for(auto &i: gseq->_reverse_rnas[1]->_exons){
//      cout<< i->_iv.left()<<endl;
//   }
//#endif
   gseq = nullptr;

}

void GffReader::sortExonOrderInMinusStrand(){
   for(const auto &i: _g_seqs){
      for(const auto &j : i->_reverse_rnas){
         if (j->_exons.size() > 1) {
            if (j->_exons[0]->_iv > j->_exons[1]->_iv) {
               reverse(j->_exons.begin(), j->_exons.end());
            }
         }
      }
   }
}
