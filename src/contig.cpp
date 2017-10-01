
#include <algorithm>
#include <iostream>
//#include <mutex>
//#include <iterator>
#include "contig.h"
#include "read.hpp"
using namespace std;



bool readhit_2_genomicFeats(const ReadHit & rh, vector<GenomicFeature> & feats){
   /*
    * Convert read hit to vector of GenomicFeature. Return false if a read is problematic.
    */
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
            LOG(WARNING)<<"Read at reference id: "<< rh.ref_id()+1 << " and position "<< rh.left()<<" has suspicious DELETION\n";
            return false;
         }
         feats.back()._match_op._len += cig[i]._length;
         offset += cig[i]._length;
         break;
      case INS:
         if(i<1 || i+1 == cig.size() || cig[i-1]._type != MATCH || cig[i+1]._type != MATCH){
            LOG(WARNING)<<"Read at reference id: "<< rh.ref_id()+1 << " and position "<< rh.left()<<" has suspicious INSERTION\n";
            return false;
         }
         break;
      case SOFT_CLIP:
         break;
      default:
         LOG(WARNING)<<"Read at reference id: "<< rh.ref_id()+1 << " and position "<< rh.left()<<" has unknown CIGAR\n";
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
}


int GenomicFeature::len() const
{
   return _match_op._len;
}

void GenomicFeature::left(uint left)
{
      int right = _genomic_offset + _match_op._len - 1;
      _genomic_offset = left;
      _match_op._len = right -left + 1;
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
   else return 0;
}

bool GenomicFeature::overlap_in_genome(const GenomicFeature & feat, const uint left, const uint right)
{
   return ( feat.left() <= right && left <= feat.right() );
}

int GenomicFeature::overlap_len_in_genome(const GenomicFeature& feat, const uint left, const uint right)
{
   //Assume in the same chromosome or contig
   if( feat.left() <= right && left <= feat.right() ) {
      return min(feat.right(), right) - max(feat.left(), left)+1;
   }
   else return 0;
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

bool GenomicFeature::compatible_2_read(const Contig& read) const
{
   if( _match_op._code != Match_t::S_MATCH) return false;
   for(auto it = read._genomic_feats.cbegin(); it != read._genomic_feats.cend(); ++it){
      if(it->_match_op._code == Match_t::S_INTRON) return false;
      if(!this->contains(*it)) return false;
   }
   return true;
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

bool operator==(const GenomicFeature &lhs, const GenomicFeature & rhs)
{
   return ( lhs._match_op._code == rhs._match_op._code &&
         lhs._genomic_offset == rhs._genomic_offset &&
         lhs._match_op._len == rhs._match_op._len);
}

bool operator!=(const GenomicFeature &lhs, const GenomicFeature &rhs)
{
   return  !(lhs == rhs);
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
    * The genomic features has to be sorted.
    */

   for(size_t i=0; i<feats.size(); ++i){
      result.push_back(feats[i]);
      GenomicFeature & f = result.back();
      while(i + 1< feats.size() &&
         f.right() + 1 == feats[i + 1].left() &&
         f._match_op._code == feats[i + 1]._match_op._code)
      {
         f._match_op._len += feats[i + 1]._match_op._len;
         ++i;
      }
   }
}


Contig::Contig(const PairedHit& ph):
      _is_ref(false),
      _ref_id(ph.ref_id()),
      _contig_id(ph.read_id()),
      _strand(ph.strand())
{
   if (!ph.is_paired()) {
      if(ph._left_read){
         _single_read_orit = SingleOrit_t::Forward;
      }
      else{
         _single_read_orit = SingleOrit_t::Reverse;
      }
   }
   vector<GenomicFeature> g_feats;
   if(ph._left_read && ph._right_read){
      readhit_2_genomicFeats(ph.left_read_obj(), g_feats);
      readhit_2_genomicFeats(ph.right_read_obj(), g_feats);
      int gap_len = (int)ph._right_read->left() - (int)ph._left_read->right() -1 ;
      if( gap_len > 0){
         g_feats.push_back(GenomicFeature(Match_t::S_GAP, ph._left_read->right()+1, (uint)gap_len));
      } else {
         std::sort(g_feats.begin(), g_feats.end());
         g_feats = merge_genomicFeats(g_feats);
      }
   }

   else{
      if(ph._right_read){
         readhit_2_genomicFeats(ph.right_read_obj(), g_feats);
      }
      if(ph._left_read){
         readhit_2_genomicFeats(ph.left_read_obj(), g_feats);
      }
   }

   if (g_feats.empty()) {
      _ref_id = -1;
      _contig_id = -1;
      _strand = Strand_t::StrandUnknown;
   } else {
      sort(g_feats.begin(), g_feats.end());
//      uint check_right = ph.left_pos();
//      for(auto i = g_feats.cbegin(); i != g_feats.cend(); ++i){
//         check_right += i->_match_op._len;
//      }
//      assert(check_right = ph.right_pos()+1);
      _genomic_feats = move(g_feats);
      _mass = ph.collapse_mass();

   }
}

//Contig::Contig(const ExonBin& eb)
//{
//   vector<GenomicFeature> g_feats;
//   g_feats.push_back(*eb._exons_in_bin.front());
//   for(size_t i = 1; i< eb._exons_in_bin.size(); ++i){
//      if(eb._exons_in_bin[i]->left()-eb._exons_in_bin[i-1]->right() == 1) {
//         continue;
//      }
//      GenomicFeature intron(Match_t::S_INTRON,
//            eb._exons_in_bin[i-1]->right()+1,
//            eb._exons_in_bin[i]->left()-eb._exons_in_bin[i-1]->right()-1);
//      g_feats.push_back(intron);
//      g_feats.push_back(*eb._exons_in_bin[i]);
//   }
//   vector<GenomicFeature> merged_feats;
//   GenomicFeature::mergeFeatures(g_feats, merged_feats);
//   _genomic_feats = move(merged_feats);
//   _strand = Strand_t::StrandUnknown;
//   _ref_id = eb.ref_id();
//   _mass = 0.0;
//   _is_ref = false;
//}


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

SingleOrit_t Contig::single_read_orit() const
{
   return _single_read_orit;
}

//bool Contig::operator<(const Contig &rhs) const
//{
//   if(_ref_id != rhs._ref_id){
//      return _ref_id < rhs._ref_id;
//   }
//   if(left() != rhs.left()){
//      return left() < rhs.left();
//   }
//   if(right() != rhs.right()){
//      return right() < rhs.right();
//   }
//   if(_genomic_feats.size() != rhs._genomic_feats.size()){
//      return _genomic_feats.size() < rhs._genomic_feats.size();
//   }
//   for(size_t i=0; i<_genomic_feats.size(); ++i){
//      if(_genomic_feats[i] != rhs._genomic_feats[i])
//         return _genomic_feats[i] < rhs._genomic_feats[i];
//   }
//   return false;
//}

bool Contig::operator<(const Contig &rhs) const {
   if(_ref_id != rhs._ref_id){
      return _ref_id < rhs._ref_id;
   }
   return _genomic_feats < rhs._genomic_feats;
}

bool operator==(const Contig &lhs, const Contig & rhs)
{
   return lhs._genomic_feats == rhs._genomic_feats;
}

bool Contig::is_single_read() const
{
   if(_is_ref) return false;
   for(auto const & gf: _genomic_feats){
      if(gf._match_op._code == S_GAP) return false;
   }
   return true;
}

uint Contig::gap_left() const
{
   if(is_single_read()) return 0;
   for(auto const & gf: _genomic_feats){
      if(gf._match_op._code == S_GAP) return gf.left();
   }
   return 0;
}

uint Contig::gap_right() const
{
   if(is_single_read()) return 0;
   for(auto const & gf: _genomic_feats){
      if(gf._match_op._code == S_GAP) return gf.right();
   }
   return 0;
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

int Contig::exonic_overlaps_len(const Contig &iso,
                                   const uint left,
                                   const uint right)
/*
 * Searching for exon fragments in a transcripts which overlap a read
 * This function should be only called if the read and the transcript are compatible
 */
{
   int collected_len = 0;
   for(auto const& gfeat : iso._genomic_feats){
      if(gfeat._match_op._code == Match_t::S_MATCH)
         collected_len +=  GenomicFeature::overlap_len_in_genome(gfeat, left, right);
   }
   return collected_len;
}

int Contig::fragment_len(const Contig& read, const Contig& iso) 
{
   if(is_compatible(read, iso))
      return exonic_overlaps_len(iso, read.left(), read.right());     
   else
      return 0;
}
   
int Contig::exonic_length() const
{
   int len = 0;
   for(auto const gf: _genomic_feats){
      if(gf._match_op._code == Match_t::S_MATCH){
         len += gf._match_op._len;
      }
   }
   return len;
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

//int Contig::infer_inner_dist(const Contig & isoform, const Contig & hit){
//   int gap_len = 0;
//   for(const GenomicFeature &hit_gf: hit._genomic_feats){
//      if(hit_gf._match_op._code == Match_t::S_GAP){
//         for(const GenomicFeature &iso_gf: isoform._genomic_feats){
//            if(iso_gf._match_op._code == Match_t::S_MATCH){
//                  gap_len += GenomicFeature::overlap_len(hit_gf, iso_gf);
//            }
//         } // inner for loop
//      }
//   }
//   return gap_len;
//}


uint Contig::read_start_from_iso(const Contig &iso, const Contig& hit)
{
   //if( !Contig::is_compatible(hit, iso) ) return 0;
   uint read_start = hit.left();
   uint dist_from_start = 0;
   for(size_t i=0; i<iso._genomic_feats.size(); ++i){
      if(iso._genomic_feats[i]._match_op._code == S_MATCH){
         if(read_start >= iso._genomic_feats[i].left() && read_start <= iso._genomic_feats[i].right() ){
            dist_from_start += read_start - iso._genomic_feats[i].left() + 1;
            break;
         }
         else if(read_start > iso._genomic_feats[i].right()){
            dist_from_start += iso._genomic_feats[i].len();
         }
         else{
            return 0;
         }
      }
   }

   return dist_from_start;
}

vector<double> Contig::start_site_dist(const Contig & iso, const vector<Contig> & hits)
{
   int nsites = iso.exonic_length();
   vector<double> start_dist(nsites, 0.0);
   vector<double> end_dist(nsites, 0.0);
   for(auto & h: hits){
      uint s = Contig::read_start_from_iso(iso, h);
      if(s == 0) continue;
      //cerr<<"nistes: " <<nsites<<endl;
      //cerr<<"s: "<<s<<endl;
      start_dist[s - 1] += h.mass();
   }
   return start_dist;
}



bool Contig::is_compatible(const Contig &read, const Contig &isoform)
/*
 * Judge if a read is compatible with a transcript.
 * Looking for first overlaping exon
 */
{
   if(read._is_ref) return false;

   vector<const GenomicFeature*> exons;
   for(size_t i = 0; i<isoform._genomic_feats.size(); ++i){
      if(isoform._genomic_feats[i]._match_op._code == S_MATCH){
         exons.push_back(&isoform._genomic_feats[i]);
      }
   }

   const GenomicFeature & first_feat = read._genomic_feats[0];
   auto first_exon = lower_bound(exons.begin(), exons.end(), first_feat,
         [](const GenomicFeature* gf, const GenomicFeature & first)
         {return gf->right() < first.left();});
   if(first_exon == exons.end()){
      return false;
   }
   else{
      if(!(*first_exon)->contains(first_feat)) return false;
      else{
         auto it = first_exon;
         for(size_t i = 1; i != read._genomic_feats.size(); ++i){
            if(read._genomic_feats[i]._match_op._code == S_GAP){
               continue;
            }
            else if(read._genomic_feats[i]._match_op._code == S_INTRON){
               size_t next_intron_offset = 2* distance(exons.begin(), it) + 1;
               if(next_intron_offset >= isoform._genomic_feats.size())
                  return false;
               if(read._genomic_feats[i] != isoform._genomic_feats[next_intron_offset])
                  return false;
            }
            else{
               assert(read._genomic_feats[i]._match_op._code == S_MATCH);
               it = find_if(it, exons.end(),
                     [&](const GenomicFeature *p)
                     {return p->contains(read._genomic_feats[i]);
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

double Contig::avg_doc() const
{
   int n = 0;
   double doc = 0.0;
   for(size_t i = 0; i<_genomic_feats.size(); ++i){
      if(_genomic_feats[i]._match_op._code == S_MATCH){
         doc += _genomic_feats[i].avg_doc();
         ++n;
      }
   }
   assert(n > 0);
   return doc/n;
}

bool Contig::is_compatible(const Contig &isoform, const GenomicFeature &feat)
{
   if(feat._match_op._code != Match_t::S_MATCH) return false;
   vector<const GenomicFeature*> exons;
   for(size_t i = 0; i<isoform._genomic_feats.size(); ++i){
      if(isoform._genomic_feats[i]._match_op._code == S_MATCH){
         exons.push_back(&isoform._genomic_feats[i]);
      }
   }
   auto first_exon = lower_bound(exons.begin(), exons.end(), feat,
      [](const GenomicFeature* gf, const GenomicFeature & first)
      {return gf->right() < first.left();});
   if(first_exon == exons.end()){
      return false;
   }
   else{
      if(!(*first_exon)->contains(feat)) return false;
   }
   return true;
}

void Contig::print2gtf(FILE *pFile,
                       const RefSeqTable &ref_lookup,
                       const string fpkm,
                       const string tpm,
                       string gene_id, string tscp_id) const {

   const char* ref = ref_lookup.ref_real_name(_ref_id).c_str();

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
   char fpkm_c[12];
   char tpm_c[12];
   strncpy(fpkm_c, fpkm.c_str(), sizeof(fpkm_c));
   fpkm_c[sizeof(fpkm_c) - 1] = 0;
   strncpy(tpm_c, tpm.c_str(), sizeof(tpm_c));
   tpm_c[sizeof(tpm_c) - 1] = 0;

   strcpy(gff_attr, "gene_id \"");
   strcat(gff_attr, gene_id.c_str());
   strcat(gff_attr, "\"; transcript_id \"");
   strcat(gff_attr, tscp_id.c_str());
   strcat(gff_attr, "\";");
   strcat(gff_attr, "FPKM \"");
   strcat(gff_attr, fpkm_c);
   strcat(gff_attr, "\";");
   strcat(gff_attr, "Frac \"");
   strcat(gff_attr, tpm_c);
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
