/*
 * contig.cpp
 *
 *  Created on: Jan 19, 2015
 *      Author: ruolin
 */



#include <algorithm>
#include <iostream>
//#include <iterator>
#include "contig.h"
#include "read.h"
using namespace std;



bool readhit_2_genomicFeats(const ReadHit & rh, vector<GenomicFeature> & feats){
   vector<CigarOp> cig = rh.cigar();
   uint offset = rh.left();
   assert(!cig.empty());
   for(size_t i = 0; i< cig.size(); ++i){
      switch(cig[i]._type)
      {
      case MATCH:
         feats.push_back(GenomicFeature(Match_t::S_MATCH, offset, cig[i]._length));
         offset += cig[i]._length;
         break;
      case REF_SKIP:
         feats.push_back(GenomicFeature(Match_t::S_INTRON, offset, cig[i]._length));
         offset += cig[i]._length;
         break;
      case DEL:
         if(i<1 || i+1 == cig.size() || cig[i-1]._type != MATCH || cig[i+1]._type != MATCH){
            LOG_INPUT("Read at reference id: ", rh.ref_id()+1, " and position ", rh.left(), " has suspicious DELETION");
            return false;
         }
         feats.back()._match_op._len += cig[i]._length;
         offset += cig[i]._length;
         ++i;
         feats.back()._match_op._len += cig[i]._length;
         offset += cig[i]._length;
         break;
      case INS:
         if(i<1 || i+1 == cig.size() || cig[i-1]._type != MATCH || cig[i+1]._type != MATCH){
            LOG_INPUT("Read at reference id: ", rh.ref_id()+1, " and position ", rh.left(), " has suspicious INSERTION");
            return false;
         }
         ++i;
         feats.back()._match_op._len += cig[i]._length;
         offset += cig[i]._length;
         break;
      case SOFT_CLIP:
         break;
      default:
         LOG_INPUT("Read at reference id: ", rh.ref_id()+1, " and position ", rh.left(), " has unknown cigar");
         return false;
      }
   }
   return true;
}

GenomicFeature::GenomicFeature(const Match_t& code, uint offset, int len):
         _genomic_offset(offset),
         _match_op (MatchOp(code, (uint32_t)len))
{
   _avg_cov=0.0;
//   if(len <= 0)
//   cout<<len<<endl;
   if(len ==0 && code != Match_t::S_GAP){
      assert (false);
   }
   assert(len >= 0);
}

void GenomicFeature::left(uint left)
{
      int right = _genomic_offset + _match_op._len;
      _genomic_offset = left;
      _match_op._len = right -left;
}

uint GenomicFeature::left() const { return _genomic_offset;}

uint GenomicFeature::right() const
{
   return _genomic_offset + _match_op._len-1;
}

void GenomicFeature::right(uint right)
{
      _match_op._len = right - _genomic_offset + 1;
}

void GenomicFeature::avg_doc(double coverage)
{
   _avg_cov = coverage;
}

double GenomicFeature::avg_doc() const
{
   return _avg_cov;
}

bool GenomicFeature::overlaps(const GenomicFeature& lhs, const GenomicFeature& rhs)
{
   //Assume in the same chromosome or contig
   return lhs.left() <= rhs.right() && rhs.left() <= lhs.right();
}

int GenomicFeature::overlap_len(const GenomicFeature &lhs, const GenomicFeature &rhs)
{
   //Assume in the same chromosome or contig
   if(GenomicFeature::overlaps(lhs, rhs)){
      return min(lhs.right(), rhs.right()) - max(lhs.left(), rhs.left());
   }
   else{
      return 0;
   }
}

bool GenomicFeature::overlap_in_genome(const GenomicFeature& feat, uint s, uint e){
   //Assume in the same chromosome or contig
   if(feat.left() >= s && feat.left() < e)
      return true;
   if(feat.right() > s && feat.right() < e)
      return true;
   if(s >= feat.left() && s < feat.right())
      return true;
   if (e > feat.left() && e < feat.right())
      return true;
   return false;
}

bool GenomicFeature::contains(const GenomicFeature& other, int small_extent) const
{
   if (left() -small_extent <= other.left() && right() + small_extent >= other.right())
      return true;
   return false;
}

bool GenomicFeature::properly_contains(const GenomicFeature& other) const
{
   if( (left() < other.left() && right() >= other.right() ) ||
         (left() <= other.left() && right() > other.right()) )
      return true;
   return false;
}

int match_length(const GenomicFeature &op, int left, int right)
{
   int len = 0;
   int left_off = op.left();
   if(left_off + op._match_op._len > left && left_off < right){
      if(left_off > left){
         if(left_off + op._match_op._len <= right + 1)
            len += op._match_op._len;
         else
            len += right -left_off;
      }
      else{
         if(left_off + op._match_op._len <= right +1)
            len += (left_off + op._match_op._len - left);
         else
            return right - left;
      }
   }
   return len;
}

bool GenomicFeature::operator==(const GenomicFeature & rhs) const
{
   return ( _match_op._code == rhs._match_op._code &&
         _genomic_offset == rhs._genomic_offset &&
         _match_op._len == rhs._match_op._len);
}

bool GenomicFeature::operator<(const GenomicFeature & rhs) const
{
   if(_genomic_offset != rhs._genomic_offset)
      return _genomic_offset < rhs._genomic_offset;
   if(_match_op._len != rhs._match_op._len)
      return _match_op._len < rhs._match_op._len;
   return false;
}


void GenomicFeature::mergeFeatures(const vector<GenomicFeature> & feats, vector<GenomicFeature> &result){
   /*
    * Merge the introduced exon segments back together
    */

   for(size_t i=0; i<feats.size(); ++i){
      result.push_back(feats[i]);
      GenomicFeature & f = result.back();
      while(i<feats.size()-1 &&
         f.right() + 1 == feats[i+1].left() &&
         f._match_op._code == feats[i+1]._match_op._code)
      {
         f._match_op._len += feats[i+1]._match_op._len;
         ++i;
      }
   }
}



Contig::Contig(RefID ref_id, Strand_t strand, const vector<GenomicFeature> &feats, double mass, bool is_ref):
      _genomic_feats(feats),
      _is_ref(is_ref),
      _ref_id(ref_id),
      _strand(strand),
      _mass(mass)
   {
      assert(_genomic_feats.front()._match_op._code == Match_t::S_MATCH);
      assert(_genomic_feats.back()._match_op._code == Match_t::S_MATCH);
   }

Contig::Contig(const PairedHit& ph):
      _is_ref(false),
      _ref_id(ph.ref_id()),
      _strand(ph.strand())
{
   vector<GenomicFeature> g_feats;
   if(ph._left_read){
      if(readhit_2_genomicFeats(ph.left_read_obj(), g_feats)){
         if(ph._right_read){
            int gap_len = (int)ph._right_read->left() - (int)ph._left_read->right() -1 ;
//            if(gap_len == 0){
//               cout<<ph.left_read_obj().left()<<"-"<<ph.right_read_obj().left()<<endl;
//               exit(0);
//            }
            if( gap_len < 0){
               gap_len = 0;
               LOG_ERR("Read at reference id: ", ph.ref_id()+1, " and position ", ph.left_pos(), " has suspicious GAP");
            }
            g_feats.push_back(GenomicFeature(Match_t::S_GAP, ph._left_read->right()+1, (uint)gap_len));
            readhit_2_genomicFeats(ph.right_read_obj(), g_feats);
         }
      }

   }

   else{
      if(ph._right_read){
         readhit_2_genomicFeats(ph.right_read_obj(), g_feats);
      }
      assert(false);
   }

   if(!is_sorted(g_feats.begin(), g_feats.end()))
   {
      sort(g_feats.begin(), g_feats.end());
   }

   uint check_right = ph.left_pos();
   for(auto i = g_feats.cbegin(); i != g_feats.cend(); ++i){
      check_right += i->_match_op._len;
   }
   assert(check_right = ph.right_pos()+1);
   _genomic_feats = move(g_feats);
   _mass = ph.collapse_mass();
}

Contig::Contig(const ExonBin& eb)
{
   vector<GenomicFeature> g_feats;
   g_feats.push_back(*eb._exons_in_bin.front());
   for(size_t i = 1; i< eb._exons_in_bin.size(); ++i){
      if(eb._exons_in_bin[i]->left()-eb._exons_in_bin[i-1]->right() == 1) {
         continue;
      }
      GenomicFeature intron(Match_t::S_INTRON,
            eb._exons_in_bin[i-1]->right()+1,
            eb._exons_in_bin[i]->left()-eb._exons_in_bin[i-1]->right()-1);
      g_feats.push_back(intron);
      g_feats.push_back(*eb._exons_in_bin[i]);
   }
   vector<GenomicFeature> merged_feats;
   GenomicFeature::mergeFeatures(g_feats, merged_feats);
   _genomic_feats = move(merged_feats);
   _strand = Strand_t::StrandUnknown;
   _ref_id = eb.ref_id();
   _mass = 0.0;
   _is_ref = false;
}


uint Contig::left() const
{
   //cout<<_genomic_feats.front().left()<<endl;

   return _genomic_feats.front().left();
}

uint Contig::right() const
{

   return _genomic_feats.back().right();
}

float Contig::mass() const
{
   return _mass;
}

void Contig::mass(float m)
{
   _mass = m;
}

bool Contig::operator<(const Contig &rhs) const
{
   if(_ref_id != rhs._ref_id){
      return _ref_id < rhs._ref_id;
   }
   if(left() != rhs.left()){
      return left() < rhs.left();
   }
   if(right() != rhs.right()){
      return right() < rhs.right();
   }
   if(_genomic_feats.size() != rhs._genomic_feats.size()){
      return _genomic_feats.size() < rhs._genomic_feats.size();
   }
   for(size_t i=0; i<_genomic_feats.size(); ++i){
      if(_genomic_feats[i] != rhs._genomic_feats[i])
         return _genomic_feats[i] < rhs._genomic_feats[i];
   }
   return false;
}

RefID Contig::ref_id() const
{
   return _ref_id;
}

Strand_t Contig::strand() const
{
   return _strand;
}

const string Contig::annotated_trans_id() const{
   return _annotated_trans_id;
}

void Contig::annotated_trans_id(string str){
   assert(!str.empty());
   _annotated_trans_id = str;
}

size_t Contig::featSize() const{
   return _genomic_feats.size();
}

bool Contig::overlaps_directional(const Contig &lhs, const Contig &rhs){
   if(lhs.ref_id() != rhs.ref_id())
      return false;
   if(lhs.strand() != rhs.strand())
      return false;
   return overlaps_locally(lhs.left(), lhs.right(), rhs.left(), rhs.right());
}

bool Contig::overlaps_only_on_exons(const Contig &ct, const GenomicFeature & exon)
/*
 * Searching for exon fragments witch overlap a read fragment
 */
{
   bool result = false;
   for(auto gfeat : ct._genomic_feats){
      if(gfeat._match_op._code != Match_t::S_INTRON){
         if(GenomicFeature::overlaps(gfeat, exon))
            result = true;
      }
   }
   return result;
}

bool Contig::is_contained_in(const Contig & small, const Contig & large)
/*
 * Judge if a exon bin is compatible with a transcript
 */
{
   // single exon case;
   if(small._genomic_feats.size() == 1){
      assert(small._genomic_feats[0]._match_op._code == Match_t::S_MATCH);
      for(size_t i = 0; i<large._genomic_feats.size(); ++i){
         if(large._genomic_feats[i]._match_op._code == Match_t::S_MATCH){
            if(large._genomic_feats[i].contains(small._genomic_feats[0])) return true;
         }
      }

      return false;
   }

   // All introns are matched  -> contained.
   vector<GenomicFeature> small_introns;
   vector<GenomicFeature> large_introns;
   for(size_t i = 0; i< small._genomic_feats.size(); ++i){
      if(small._genomic_feats[i]._match_op._code == Match_t::S_INTRON){
         small_introns.push_back(small._genomic_feats[i]);
      }
   }
   for(size_t i = 0; i< large._genomic_feats.size(); ++i){
      if(large._genomic_feats[i]._match_op._code == Match_t::S_INTRON){
         large_introns.push_back(large._genomic_feats[i]);
      }
   }

   vector<bool> res(small_introns.size(), false);
   for(size_t i = 0; i< small_introns.size(); ++i){
      if(binary_search(large_introns.begin(), large_introns.end(), small_introns[i])){
         res[i] = true;
      }
   }

   for(size_t i = 0; i< res.size(); ++i){
      if(res[i] == false) return false;
   }
   return true;
}

int Contig::infer_insert_size(const Contig & isoform, const Contig & hit){
   int gap_len = 0;
   for(const GenomicFeature &hit_gf: hit._genomic_feats){
      if(hit_gf._match_op._code == Match_t::S_GAP){
         for(const GenomicFeature &iso_gf: isoform._genomic_feats){
            if(iso_gf._match_op._code == Match_t::S_MATCH){
                  gap_len += GenomicFeature::overlap_len(hit_gf, iso_gf);
            }
         } // inner for loop
      }
   }
   return gap_len;
}

bool Contig::is_compatible(const Contig &read, const Contig &isoform)
/*
 * Judge if a read is compatible with a transcript.
 * Looking for first overlaping exon
 */
{
   vector<const GenomicFeature*> exons;
   for(size_t i = 0; i<isoform._genomic_feats.size(); ++i){
      if(isoform._genomic_feats[i]._match_op._code == S_MATCH){
         exons.push_back(&isoform._genomic_feats[i]);
      }
   }

   const GenomicFeature & first_feat = read._genomic_feats[0];
   auto first = lower_bound(exons.begin(), exons.end(), first_feat,
         [](const GenomicFeature* gf, const GenomicFeature & first)
         {return gf->right() < first.left();});
   if(first == exons.end()){
      return false;
   }
   else{
      if(!(*first)->contains(first_feat,2)) return false;
      else{
         auto it = first;
         for(size_t i = 1; i != read._genomic_feats.size(); ++i){
            if(read._genomic_feats[i]._match_op._code == S_GAP){
               continue;
            }
            else if(read._genomic_feats[i]._match_op._code == S_INTRON){
               size_t next_intron_offset = 2* distance(exons.begin(), it) + 1;
               if(read._genomic_feats[i] != isoform._genomic_feats[next_intron_offset]){
                  return false;
               }
            }
            else{
               assert(read._genomic_feats[i]._match_op._code == S_MATCH);
               it = find_if(it, exons.end(),
                     [&](const GenomicFeature *p)
                     {return p->contains(read._genomic_feats[i],2);
                     });
               if(it == exons.end()){
                  return false;
               }
               //cout<<"read: "<<read.left()<<": "<<read._genomic_feats[i].left()<<"it "<<(*it)->left()<<endl;
            }
         } // end foor loop
      }
   }
   return true;
}



void Contig::print2gtf(FILE *pFile, const RefSeqTable &ref_lookup, int gene_id, int tscp_id){
   const char* ref = ref_lookup.ref_name(_ref_id).c_str();
   char strand = 0;
   switch(_strand){
   case Strand_t::StrandPlus:
      strand = '+';
      break;
   case Strand_t::StrandMinus:
      strand = '-';
      break;
   default:
      strand = '.';
      break;
   }

   char gff_attr[200];
   char locus[13] = "gene.";
   char tscp[24] = "transcript.";
   char gene_str[7];
   char tscp_str[4]; // at most 999 transcripts for a gene.
   Sitoa(gene_id, gene_str, 10);
   Sitoa(tscp_id,tscp_str, 10);

   strcat(locus, gene_str);

   strcat(tscp, gene_str);
   strcat(tscp, ".");
   strcat(tscp, tscp_str);

   strcpy(gff_attr, "gene_id \"");
   strcat(gff_attr, locus);
   strcat(gff_attr, "\"; transcript_id \"");
   strcat(gff_attr, tscp);
   strcat(gff_attr, "\";");

   fprintf(pFile, "%s\t%s\t%s\t%d\t%d\t%d\t%c\t%c\t%s\n", \
         ref, "Strawberry", "transcript", left(), right(), 1000, strand, '.', gff_attr);

   int exon_num = 0;
   for(auto gfeat : _genomic_feats){
      if(gfeat._match_op._code == Match_t::S_MATCH){
         ++exon_num;
         char exon_gff_attr[200];
         char exon_id[5];
         Sitoa(exon_num, exon_id, 10);
         strcpy(exon_gff_attr, gff_attr);
         strcat(exon_gff_attr, " exon_id \"");
         strcat(exon_gff_attr, exon_id);
         strcat(exon_gff_attr, "\";");
         fprintf(pFile, "%s\t%s\t%s\t%d\t%d\t%d\t%c\t%c\t%s\n", \
            ref, "Strawberry", "exon", gfeat.left(), gfeat.right(), 1000, strand, '.', exon_gff_attr);
      }
   }

}

ExonBin::ExonBin(RefID ref_id)
{
   _ref_id = ref_id;
   _read_num  = 1;
}

RefID ExonBin::ref_id() const
{
   return _ref_id;
}

bool ExonBin::insert_exon(const GenomicFeature *exon)
{
   bool status;
   set<uint>::iterator it_ret = _coords.find(exon->left());
   if(it_ret == _coords.end()){
      _coords.insert(exon->left());
      _coords.insert(exon->right());
      _exons_in_bin.push_back(exon);
      status = true;
   }
   else{
      assert(_coords.find(exon->right()) != _coords.end());
      status = false;
   }
   return status;
}

void ExonBin::add_isoform(const vector<Isoform> &isoforms){

   for(size_t i = 0; i< isoforms.size(); ++i){
      bool compatible = true;
      //cout<<this->_exons_in_bin.size()<<endl;
      Contig exons(*this);

      if(Contig::is_contained_in(exons, isoforms[i]._contig)){
         _parent_isoforms.push_back(isoforms[i]._isoform_id);
      }
   }
}

uint ExonBin::left_most() const
{
   return *_coords.begin();
}


uint ExonBin::right_most() const
{
   return *_coords.rbegin();
}


void ExonBin::read_num_increase_by_1()
{
   ++_read_num;
}

int ExonBin::read_num() const
{
   return _read_num;
}

void ExonBin::add_hit(const Contig* hit)
{
   _hits.push_back(hit);
}

void ExonBin::add_hit(const ExonBin& exb)
{
   assert(exb._hits.size() == 1);
   _hits.push_back(exb._hits.front());
}


Isoform::Isoform(Contig contig, string gene, string iso_id):
      _contig(contig), _gene_id(gene), _isoform_id(iso_id)
{
   _bais_factor = 0.0;
}
