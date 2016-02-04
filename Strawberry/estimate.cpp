/*
>HEADER
    Copyright (c) 2015 Ruolin Liu rliu0606@gmail.com
    This file is part of Strawberry.
    Strawberry is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Strawberry is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Strawberry.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
*/




#include <iostream>
#include <cassert>
#include <iterator>
#include <random>
#include <Eigen/Dense>
#include <stdexcept>
#include "estimate.h"
#include "fasta.h"
#include "contig.h"

// choose exact integral type

// program and solution types

const double Estimation::_kMinTPM = kMinIsoformFrac * 1e6;
//const double Estimation::_kMinTPM = 0;

uint generate_pair_end(const Contig& ct, const Contig& orig_read, int dist_need_2_extend,  SingleOrit_t orit, Contig & inferred_mp){
/*
 * Given one end of a read pair, its orientation, simulated inner distance and its possible transcript origin
 * Genearte a simulated pair-end read as Contig inferred_mp
 *
 * Find the exon that the read inner part hit and extend the read from that position to the end of the exon
 * The extended distance is calculated as the total distance need_to_extend.
 * If still_need is negative, roll back. Stop till still_need is 0.
 *
 *       isoform:                      |@@@@@@@....@@@@@@@@|
 *       raw read, need to extend 4bp     @@@@
 * =====================================================================
 *    1. extend 1, still_need 3           @@@@@
 *    2. extend 8, still_need -5          @@@@@....@@@@@@@@
 *    3. extend -5, still_need 0          @@@@@....@@@
 */
   vector<GenomicFeature> exons;
   inferred_mp = orig_read;
   inferred_mp.mass(orig_read.mass());

   for(auto gf:ct._genomic_feats){
      if(gf._match_op._code == Match_t::S_MATCH){
         exons.push_back(gf);
      }
   }

   if(orit == SingleOrit_t::Forward){
      uint start = inferred_mp.right();

      auto it = lower_bound(exons.cbegin(), exons.cend(), start,
         [](const GenomicFeature &gf, uint s){return gf.right() < s;});
      assert(it != exons.cend());
      dist_need_2_extend = dist_need_2_extend + start  - it->right();
      while(dist_need_2_extend >0 && ++it != exons.cend()){
         if(inferred_mp._genomic_feats.back()._match_op._code == Match_t::S_INTRON){
             GenomicFeature gf(Match_t::S_MATCH, start, (it-1)->right() -start +1);
             inferred_mp._genomic_feats.push_back(gf);
         }
         else{
            inferred_mp._genomic_feats.back().right( (it-1)->right() );
         }
         GenomicFeature gf(Match_t::S_INTRON, (it-1)->right()+1, it->left() - 1 - (it-1)->right());
         inferred_mp._genomic_feats.push_back(gf);
         start = it->left();
         dist_need_2_extend = dist_need_2_extend + start - it->right();
      }
      if(it == exons.cend() && dist_need_2_extend >0)
         return 0;

      if(inferred_mp._genomic_feats.back()._match_op._code == Match_t::S_INTRON){
         GenomicFeature gf(Match_t::S_MATCH, it->left(), it->right() + dist_need_2_extend - it->left() +1);
         inferred_mp._genomic_feats.push_back(gf);
      }
      else{
         inferred_mp._genomic_feats.back().right(dist_need_2_extend + it->right());
      }

      return dist_need_2_extend + it->right();
   }


   if(orit == SingleOrit_t::Reverse){
      uint start = inferred_mp.left();
      auto it = lower_bound(exons.crbegin(), exons.crend(), start,
            [](const GenomicFeature &gf, uint s){return gf.left() > s;});
      assert(it != exons.crend());
      dist_need_2_extend = dist_need_2_extend - start + it->left();
      while(dist_need_2_extend > 0 && ++it != exons.crend()){
         if(inferred_mp._genomic_feats.front()._match_op._code == Match_t::S_INTRON){
            GenomicFeature gf(Match_t::S_MATCH, (it-1)->left(), start - (it-1)->left() +1);
            inferred_mp._genomic_feats.insert(inferred_mp._genomic_feats.begin(), gf);
         }
         else{
            inferred_mp._genomic_feats.front().left( (it-1)->left() );
         }
         GenomicFeature gf(Match_t::S_INTRON, it->right()+1, (it-1)->left() -1 - it->right());
         inferred_mp._genomic_feats.insert(inferred_mp._genomic_feats.begin(), gf);
         start = it->right();
         dist_need_2_extend = dist_need_2_extend - start + it->left();
      }
      if(it == exons.crend() && dist_need_2_extend > 0)
         return 0;

      if(inferred_mp._genomic_feats.front()._match_op._code == Match_t::S_INTRON){
         GenomicFeature gf(Match_t::S_MATCH, it->left()-dist_need_2_extend, it->right() - it->left() + dist_need_2_extend +1);
         inferred_mp._genomic_feats.insert(inferred_mp._genomic_feats.begin(), gf);
      }
      else{
         inferred_mp._genomic_feats.front().left(it->left() - dist_need_2_extend);
      }
      return it->left() - dist_need_2_extend;
   }
   return 0;
}


int no_gap_ef(const int l_left, const int l_right, const int l_int, const int fl)
/*
 * Effective length of a fragment hitting both ends but
 * do not consider whether it hits the inner segments.
 */
{
   if(fl < l_int + 2) return 0;
   if(fl > l_left + l_right + l_int) return 0;
   int mid = fl - l_int -1;
   return min(l_left, mid) + min(l_right, mid) - mid;
}

int gap_ef(const int l_left, const int l_right, const int l_int, const int rl, const int gap)
/*
 * Effective length of a fragment hitting both ends but
 * CANNOT hit any inner segmetns!
 * This time we need to consider the read gap between two ends.
 */

{
   if(2*rl+gap < l_int +2) return 0;
   if(2*rl+gap > l_left + l_right + l_int) return 0;
   int start = max(rl, l_left + l_int - gap -1);
   int end = min(l_left, l_left+l_right+l_int - gap - rl);
   return ( max(0, end-start));
}


ExonBin::ExonBin(std::set<pair<uint,uint>> coordinate): _coords(coordinate){}

bool ExonBin::operator==(const ExonBin& rhs) const
{
   return _coords == rhs._coords;
}

//RefID ExonBin::ref_id() const
//{
//   return _ref_id;
//}
//
//bool ExonBin::insert_exon(const GenomicFeature *exon)
//{
//   bool status;
//   set<uint>::iterator it_ret = _coords.find(exon->left());
//   if(it_ret == _coords.end()){
//      _coords.insert(exon->left());
//      _coords.insert(exon->right());
//      _exons_in_bin.push_back(exon);
//      status = true;
//   }
//   else{
//      assert(_coords.find(exon->right()) != _coords.end());
//      status = false;
//   }
//   return status;
//}
//
//void ExonBin::add_isoform(const vector<Isoform> &isoforms){
//
//   for(size_t i = 0; i< isoforms.size(); ++i){
//      bool compatible = true;
//      //cout<<this->_exons_in_bin.size()<<endl;
//      Contig exons(*this);
//
//      if(Contig::is_contained_in(exons, isoforms[i]._contig)){
//         _parent_isoforms.push_back(isoforms[i]._isoform_id);
//      }
//   }
//}

uint ExonBin::left() const
{
   return _coords.cbegin()->first;
}


uint ExonBin::right() const
{
   return _coords.crbegin()->second;
}


void ExonBin::add_read_mass(float mass)
{
   _whole_read_mass += mass;
}

bool ExonBin::add_frag(const Contig& fg)
{
   auto ret = _frags.insert(fg);
#ifdef DEBUG
//   if(ret.second == false){
//      cout<<(*ret.first).left()<<" has existed"<<endl;
//      cout<<fg.left()<<" want to be added"<<endl;
//   }
#endif
   return true;
}

int ExonBin::num_exons() const
{
   return _coords.size();
}

float ExonBin::read_count() const
{
   float sum = 0.0 ;
   for(auto it = _frags.cbegin(); it != _frags.cend(); ++it ){
//#ifdef DEBUG
//      cout<<"it->left "<<it->left()<<" mass: "<<it->mass()<<endl;
//#endif
      sum += it->mass();
   }
   return sum;
}

int ExonBin::bin_len() const
{
   int bin_len = 0;
   for(auto c = _coords.cbegin(); c != _coords.cend(); ++c){
      bin_len = c->second - c->first;
   }
   return bin_len;
}

double ExonBin::GC_content() const
{

}

double ExonBin::avg_frag_len() const
{
   assert(_frags.empty() == false);
   int sum = 0;
   for(auto it = _frags.cbegin(); it != _frags.cend(); ++it){
      sum += it->exonic_length();
   }
   return (double) sum / _frags.size();
}

RefID ExonBin::ref_id() const
{
   assert(_frags.empty() == false);
   return _frags.cbegin()->ref_id();
}


vector<uint> ExonBin::bin_under_iso(const Isoform& iso,
         vector<pair<uint, uint>> & exon_coords) const
/*
 * NOTE: The number of segments of a bin might be less than the number of
 * segments under its compatible isoform. This is because the fragment gap.
 *
 * Return two things.
 * First, the implicit whole exon segments under a isoform
 * Second, the idx of the exon segments that are implicit.
 */
{
   vector<uint> idx;
   const vector<GenomicFeature> & exons = iso._exon_segs;
   vector<uint> start_pos(exons.size());
   for(uint i=0; i<exons.size(); ++i){
      start_pos[i] = exons[i].left();
   }
   vector<uint>::const_iterator low = lower_bound(start_pos.cbegin(), start_pos.cend(), this->left());

   assert(low != start_pos.cend());
   vector<uint>::const_iterator up  = lower_bound(start_pos.cbegin(), start_pos.cend(), _coords.crbegin()->first );
   assert(up  != start_pos.cend());
   for(vector<uint>::const_iterator it = low; it != up; ++it){
      const GenomicFeature & cur_exon = exons[distance(start_pos.cbegin(),it)];
      exon_coords.push_back( pair<uint,uint>( cur_exon.left(), cur_exon.right()));
   }
   const GenomicFeature & cur_exon = exons[distance(start_pos.cbegin(),up)];
   exon_coords.push_back(  pair<uint,uint>( cur_exon.left(), cur_exon.right()));

   auto c = _coords.cbegin();
   ++c;
   uint i = 1;

   while(i< exon_coords.size()-1){
      if(exon_coords[i].first < c->first){
         idx.push_back(i);
         ++i;
      }
      else if(exon_coords[i].first == c->first){
         ++i;
         ++c;
      }
      else{

         cout<<"false"<<endl;
         for(auto e: exons)
            cout<<e.left()<<"-"<<e.right()<<endl;
         cout<<"coordiantes"<<endl;
         assert(false);
      }
   }
   return idx;
}

int ExonBin::left_exon_len() const
{
   return _coords.cbegin()->second - _coords.cbegin()->first + 1;
}






int ExonBin::effective_len(const vector<uint> & seg_lens,
         const vector<uint>& implicit_idx,
         const int fl,
         const int rl
         ) const
{

   int gap = fl - 2*rl;

   if(seg_lens.size() == 1){ // single exon
      return seg_lens[0] - fl + 1;
   }
   else if(seg_lens.size() == 2){ //double seg_lens
      return no_gap_ef(seg_lens[0], seg_lens[1], 0, fl);
   }


   else if(seg_lens.size() == 3){ // triple seg_lens
      if(implicit_idx.size() == 1){
         return gap_ef(seg_lens[0], seg_lens[2], seg_lens[1], rl, gap);
      }

      else if(implicit_idx.size() == 0){
         return ( no_gap_ef(seg_lens[0], seg_lens[2], seg_lens[1], fl) - gap_ef(seg_lens[0], seg_lens[2], seg_lens[1], rl, gap));
      }

      else{
         assert(false);
      }
   }
   else if(seg_lens.size() == 4){ // quadruple seg_lens
      int hit14 = gap_ef(seg_lens[0], seg_lens[3], seg_lens[2] + seg_lens[1], rl, gap);
      int hit24 = gap_ef(seg_lens[3], seg_lens[1], seg_lens[2], rl ,gap);
      int hit124 = gap_ef(seg_lens[0]+seg_lens[1], seg_lens[3], seg_lens[2], rl, gap);
      int hit13 = gap_ef(seg_lens[0], seg_lens[2], seg_lens[1], rl, gap);
      int hit134 = gap_ef(seg_lens[0], seg_lens[2]+seg_lens[3], seg_lens[1], rl, gap);

      if(implicit_idx.size() == 0){
         int hit_all_124 = hit124 - hit14 - hit24;
         int hit_all_134 = hit134 - hit14 - hit13;
         int total = no_gap_ef(seg_lens[0], seg_lens[3], seg_lens[1]+seg_lens[2], fl);
         return (total - hit_all_124 - hit_all_134 - hit14);
      }
      else if(implicit_idx.size() == 2){
         return hit14;
      }
      else{
         if(implicit_idx[0] == 1){
            int hit_all_134 = hit134 - hit14 - hit13;
            return hit_all_134;
         }
         else{
            int hit_all_124 = hit124 - hit14 - hit24;
            return hit_all_124;
         }
      }
   }
   else{ // more than quadruple seg_lens
      uint num_inners = seg_lens.size() - 2;
      uint num_pos = 0;
      uint target = pow(2.0, seg_lens.size()) -1;
      for(uint idx : implicit_idx){
         target &= ~(1u << idx);
      }
      for( int i = 1; i != seg_lens[0]+1; ++i){
         uint hit = 1;
         int bp_last = fl - i - accumulate(seg_lens.begin()+1, seg_lens.end()-1,0);
         if(bp_last > *seg_lens.rbegin()) continue;
         if(bp_last < 0) assert(false);
         if(bp_last == 0) break;

         hit |= (1u << (seg_lens.size() -1));//
         //right-end cover
         int last_rest_bp = rl - bp_last;
         uint j = num_inners;
         while(last_rest_bp > 0 && j > 0){
            hit |= (1u << j);
            last_rest_bp = last_rest_bp - seg_lens[j];
            j = j-1;
         }

         //left-end cover
         int first_rest_bp = rl - i;
         j = 1;
         while(first_rest_bp > 0 && j <= num_inners){
            hit |= (1u << j);
            first_rest_bp = first_rest_bp - seg_lens[j];
            j = j+1;
         }

         if(hit == target) num_pos++;
      }
//#ifdef DEBUG
//      cout<<"fl: "<<fl<<"  num pos: "<<num_pos<<endl;
//#endif
      return num_pos;
   }
   return 0;
}
//double ExonBin::read_depth() const
/*
 * This function contains bugs.
 */
//{
//   int bin_len = this->bin_len();
//   double bin_cov = 0;
//   for(auto const & kv:_iso_2_frag_lens){
//      for(auto const &frag: kv.second){
//         bin_cov += frag.first * frag.second;
//      }
//   }
//
//   return bin_cov/bin_len;
//}


void ExonBin::add_frag_len(const int iso, const int frag_len, const float mass)
{
   add_read_mass(mass);
   auto ret = _iso_2_frag_lens.find(iso);
   if(ret == _iso_2_frag_lens.end()){
      vector<pair<int,float>> vec_lens = { pair<int,float>(frag_len, mass) };
      _iso_2_frag_lens.emplace(iso, vec_lens);
   }
   else{
      ret->second.push_back(pair<int,float>(frag_len, mass));
   }
}


Isoform::Isoform(const vector<GenomicFeature>& exons, Contig contig, int gene, int iso_id):
      _contig(contig), _gene_id(gene), _isoform_id(iso_id)
{
   for(uint i = 0; i< exons.size(); ++i){
      if(Contig::is_compatible(_contig, exons[i])){
         _exon_segs.push_back(exons[i]);
      }
   }
   _bais_factor = 0.0;
}

Estimation::Estimation(shared_ptr<InsertSize> insert_size,
                       int read_len,
                       FILE* tracker):
   _insert_size_dist(insert_size), _read_len(read_len), _p_log_file(tracker) {}

void Estimation::set_maps(
             const int& iso_id,
             const int& fg_len,
             const float& mass,
             const Contig& read,
             const set<pair<uint, uint>>& coords,
             map<set<pair<uint,uint>>, ExonBin> & exon_bin_map,
             map<int, set<set<pair<uint,uint>>>> &iso_2_bins_map)
{
   if(coords.empty()) return;
   auto find_it = exon_bin_map.find(coords);
   if(find_it == exon_bin_map.end()){
      ExonBin eb(coords);
      auto ret = exon_bin_map.emplace(coords, eb);
      ret.first->second.add_frag(read);
      ret.first->second.add_frag_len(iso_id, fg_len, mass);
   }
   else{
      find_it->second.add_frag_len(iso_id, fg_len, mass);
      find_it->second.add_frag(read);
   }
   auto find_it2 = iso_2_bins_map.find(iso_id);
   if (find_it2 == iso_2_bins_map.end()){
      set<set<pair<uint, uint>>> exon_bins = {coords};
      iso_2_bins_map.emplace(iso_id, exon_bins);
   }
   else{
      find_it2->second.insert(coords);
   }
}


void Estimation::overlap_exons(const vector<GenomicFeature>& exons,
                     const Contig& read,
                     set<pair<uint,uint>> &coords)
{
   for(auto const& gfeat : exons){
      if (gfeat._match_op._code != Match_t::S_MATCH) continue;

      for(auto const &read_f: read._genomic_feats){
         if (read_f._match_op._code != Match_t::S_MATCH) continue;

         if (GenomicFeature::overlaps(read_f, gfeat)){
            auto ret = coords.insert(pair<uint,uint>(gfeat.left(),gfeat.right()));
            //assert(ret.second)
         }
      }
   }
}



void Estimation::assign_exon_bin(
/*
 * assign reads and transcripts to exon bin.
 */
      const vector<Contig> &hits,
      const vector<Isoform> &transcripts,
      const vector<GenomicFeature> & exon_segs,
      map<set<pair<uint,uint>>, ExonBin> & exon_bin_map,
      map<int, set<set<pair<uint,uint>>>> &iso_2_bins_map)
{
   for(auto mp = hits.cbegin(); mp != hits.cend(); ++mp){

      map<set<uint>, pair<set<int>, int>> frag_mult_exonbin;
      double sr_fg_len = 0.0;
      if(mp->is_single_read() && infer_the_other_end){
         //random_device rd;
         //mt19937 gen(rd());
         mt19937 gen(3); // we use a fixed seed to make sure the output are the same every time.
         double mean = _insert_size_dist->_mean;
         double sd = _insert_size_dist->_sd;
         normal_distribution<> nd(mean, sd);

         while( (sr_fg_len = nd(gen)) <= 0){}
      }

      for(auto iso = transcripts.cbegin(); iso != transcripts.cend(); ++iso){
         if(Contig::is_compatible(*mp, iso->_contig)){
            set<pair<uint,uint>> coords;
            int frag_len = 0;

            /*For singleton, we random generate the other end */
            if(mp->is_single_read()){
               if(infer_the_other_end){
                  Contig new_read;
                  if(mp->single_read_orit() == SingleOrit_t::Forward){
                     uint other_end = generate_pair_end(iso->_contig, *mp, (int)sr_fg_len - _read_len, SingleOrit_t::Forward, new_read);
                     if(other_end == 0) continue;
                     overlap_exons(exon_segs, new_read, coords);
//#ifdef DEBUG
//   cout<<"Forward: "<<new_read.left()<<":"<<other_end<<endl;
//   cout<<"Forward: "<<new_read.right()<<endl;
//   cout<<"iso "<<iso->_isoform_id<<endl;
//   cout<<"coords";
//   for(auto const& c : coords)
//      cout<<" "<<c;
//   cout<<endl;
//#endif
                  }
                  else if (mp->single_read_orit() == SingleOrit_t::Reverse) {

                     uint other_end = generate_pair_end(iso->_contig, *mp, (int)sr_fg_len - _read_len, SingleOrit_t::Reverse, new_read);
                     if(other_end == 0) continue;
                     assert(coords.empty());
                     overlap_exons(exon_segs, new_read, coords);
//#ifdef DEBUG
//   cout<<"raw read right "<<mp->right()<<endl;
//   cout<<"Reverse: "<<other_end<<":"<<new_read.right()<<endl;
//   for(auto gf: new_read._genomic_feats){
//      cout<<gf._match_op._code<<":"<<gf.left()<<"-"<<gf.right()<<endl;
//   }
//   cout<<"Reverse: "<<new_read.left()<<endl;
//   cout<<"iso "<<iso->_isoform_id<<endl;
//   cout<<"coords";
//   for(auto const& c : coords)
//      cout<<" "<<c;
//   cout<<endl;
//#endif

                  }
                  frag_len = Contig::exonic_overlaps_len(iso->_contig, new_read.left(), new_read.right());
                  set_maps(iso->_isoform_id, frag_len, new_read.mass(), new_read, coords, exon_bin_map , iso_2_bins_map);
               } // end if infer the other
               else{
                  overlap_exons(exon_segs, *mp, coords);
                  frag_len = Contig::exonic_overlaps_len(iso->_contig, mp->left(), mp->right());
                  set_maps(iso->_isoform_id, frag_len, mp->mass(), *mp, coords, exon_bin_map , iso_2_bins_map);
               }
            } // and and single end

            else{
               overlap_exons(exon_segs,*mp, coords);
               frag_len = Contig::exonic_overlaps_len(iso->_contig, mp->left(), mp->right());
               set_maps(iso->_isoform_id, frag_len, mp->mass(), *mp, coords, exon_bin_map , iso_2_bins_map);
            }

            //assert(!coords.empty());
//#ifdef DEBUG
//            cout<<"iso "<<iso->_isoform_id<<endl;
//            cout<<"coords";
//            for(auto const& c : coords)
//               cout<<" "<<c.first<<"-"<<c.second;
//            cout<<" hit: "<<mp->left()<<" mass: "<<mp->mass()<<endl;
//            //cout<<frag_len<<""<<endl;
//#endif
         } // end if compatible condition
      }// end second inner for loop
   }
}

void Estimation::calculate_raw_iso_counts(const map<int, set<set<uint>>> &iso_2_bins_map,
          const map<set<uint>, ExonBin> & exon_bin_map)
{

}

void Estimation::calculate_bin_bias( map<set<pair<uint,uint>>, ExonBin> & exon_bin_map,
                                     shared_ptr<FaSeqGetter> &fa_getter){
   for(auto it = exon_bin_map.begin(); it != exon_bin_map.end(); ++it){
      ExonBin& eb = it->second;
      cout<<eb.ref_id()<<endl;
      cout<<eb.avg_frag_len()<<endl;
      cout<<fa_getter->fetchSeq(eb.left(),eb.right()-eb.left()+1)<<endl;
   }
}


void Estimation::theory_bin_weight(const map<int, set<set<pair<uint,uint>>>> &iso_2_bins_map,
                          const map<int,int> &iso_2_len_map,
                          const vector<Isoform>& isoforms,
                          map<set<pair<uint,uint>>, ExonBin> & exon_bin_map
                          )
{

   for(auto it = iso_2_bins_map.cbegin(); it != iso_2_bins_map.cend(); ++it){
      for(auto const &bin_coord: it->second){
         double weight = 0.0;
         assert(isoforms[it->first-1]._isoform_id == it->first);
         vector<pair<uint,uint>> exon_segs;
         vector<uint> implicit_exon_idx = exon_bin_map.at(bin_coord).bin_under_iso(isoforms[it->first-1], exon_segs);
//#ifdef DEBUG
//         cout<<"exon segments under isoform "<<it->first<<endl;
//         for(auto e: exon_segs){
//            cout<<e.first<<"-"<<e.second<<endl;
//         }
//         cout<<"in which, exon id ";
//         for(auto id: implicit_exon_idx){
//            cout<<id<<",";
//         }
//         cout<<"are implicited."<<endl;
//#endif
         vector<uint> seg_lens(exon_segs.size());
         for(uint i=0; i< seg_lens.size(); ++i){
            seg_lens[i] = exon_segs[i].second - exon_segs[i].first + 1;
         }

         int lmax = accumulate(seg_lens.begin(), seg_lens.end(),0);
         int lmin;
         if(_insert_size_dist->_use_emp)
            lmin = _insert_size_dist->_start_offset;
         else
            lmin = 2*_read_len;
         if(seg_lens.size() > 2)
            lmin = max(lmin, accumulate(seg_lens.begin() +1, seg_lens.end()-1, 0));

         for(int fl = lmin ; fl <= lmax; ++fl){
            double le_eff = exon_bin_map.at(bin_coord).effective_len(seg_lens,implicit_exon_idx, fl, _read_len);
            double tmp = _insert_size_dist->emp_dist_pdf(fl)* le_eff / (iso_2_len_map.at(it->first)-fl+1);
            weight += tmp;
         }
         exon_bin_map.at(bin_coord)._bin_weight_map[it->first] = weight;

//#ifdef DEBUG
//   cout<<"iso "<<it->first<<endl;
//   cout<<"coords";
//   for(auto const& c : bin_coord)
//      cout<<" "<<c.first<<"-"<<c.second;
//   cout<<endl;
//         cout<<"weight "<<weight<< " lmax "<<lmax<<endl;
//#endif
      }
   }

//#ifdef DEBUG
//     for(auto & il: iso_2_len_map){
//            cout<<"isoform length for "<<il.first<<": "<<il.second<<endl;
//   }
//#endif
}


void Estimation::empirical_bin_weight( const map<int, set<set<pair<uint,uint>>>> &iso_2_bins_map,
                                  const map<int,int> &iso_2_len_map,
                                  const int m,
                                  map<set<pair<uint,uint>>, ExonBin> & exon_bin_map)
{
   map<int, double> iso_weight_factor_map;
   for(auto it = iso_2_bins_map.cbegin(); it != iso_2_bins_map.cend(); ++it){
      double iso_weight = 0;
      for(auto const &bin_coord : it->second){
         double weight = 0.0;
         //iso_num = exon_bin_map.at(bin_coord)
         for(auto const &len_mass_pair: exon_bin_map.at(bin_coord)._iso_2_frag_lens.at(it->first)){
            int num_start_pos = iso_2_len_map.at(it->first)-len_mass_pair.first+1;
            weight += len_mass_pair.second * _insert_size_dist->emp_dist_pdf(len_mass_pair.first) / num_start_pos;
         }
         exon_bin_map.at(bin_coord)._bin_weight_map[it->first] = weight;
         iso_weight += weight;
      }
      iso_weight_factor_map[it->first] = iso_weight;
      cout<<"bias iso "<<iso_weight<<endl;
   }

   for(auto it = exon_bin_map.begin(); it != exon_bin_map.end(); ++it){
      for(auto iso = iso_weight_factor_map.cbegin(); iso != iso_weight_factor_map.cend(); ++iso){
         it->second._bin_weight_map[iso->first]  /= iso->second;
      }
   }

   for(auto & il: iso_2_len_map){
#ifdef DEBUG
            cout<<"isoform length for "<<il.first<<": "<<il.second<<endl;
#endif
   }

}

bool Estimation::estimate_abundances(map<set<pair<uint,uint>>, ExonBin> & exon_bin_map,
                     const double mass,
                     map<int, int>& iso_2_len_map,
                     vector<Isoform>& isoforms,
                     bool with_bias_correction,
                     shared_ptr<FaSeqGetter> &fa_getter)
{
   size_t nrow = exon_bin_map.size();
   size_t niso = isoforms.size();

   vector<int> n (nrow);
   vector<vector<double>> alpha (nrow, vector<double>(niso));
   unsigned i = 0;
   for(auto const& bin: exon_bin_map){
      n[i] = bin.second.read_count();
      for(unsigned j = 0; j< niso; ++j){
         auto ret = bin.second._bin_weight_map.find(j+1);
         if (ret == bin.second._bin_weight_map.end())
            alpha[i][j] = 0.0;
         else
            alpha[i][j] = ret->second;
      } // //inner loop
      ++i;
   } // outer loop
//#ifdef DEBUG
//   for(uint i=0; i!= alpha.size();++i){
//      for(uint j=0; j!= alpha[i].size(); ++j)
//         cout<<alpha[i][j]<<" ";
//      cout<<endl;
//   }
//#endif
   bool success;
   EmSolver em;
   if(with_bias_correction){
      calculate_bin_bias(exon_bin_map, fa_getter);
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X;
      X.setRandom(nrow, 3);
      X = X+Eigen::MatrixXd::Ones(nrow,3);
      vector<vector<double>> bias(nrow, vector<double>(3,0.0));
      for(int i=0; i< nrow; ++i)
         for(int j =0; j< 3; ++j) bias[i][j] = X(i,j);

      success = em.init(niso, n, alpha, bias);
      if(success) success = em.run();
   }
   else{
      success = em.init(niso, n, alpha);
      if(success) success = em.run();
   }
   if(success){
      for(uint i=0; i<niso; ++i){
         fprintf(_p_log_file, "isoform %d has %f raw read count.\n", i+1, em._theta[i]);
      }
      double sum_fpkm = 0.0;
      for(uint i=0; i< niso; ++i){
         uint id = isoforms[i]._isoform_id;
         double kb = 0.0;
         if(effective_len_norm){
            kb = iso_2_len_map[id] - _insert_size_dist->_mean;
            if (kb < 0){
               isoforms[i]._FPKM_s = "nan";
               continue;
            }
            kb =1e3/kb;
         }
         else{
            kb = 1e3/iso_2_len_map[id];
         }
         double rpm = 1e6/mass;
         double fpkm = em._theta[i]*rpm*kb;
         isoforms[i]._FPKM = fpkm;
         sum_fpkm += fpkm;
         isoforms[i]._FPKM_s = to_string(fpkm);
      }
      for(uint i=0; i< niso; ++i){
         if(isoforms[i]._FPKM_s == "nan"){
            isoforms[i]._TPM_s = "nan";
            continue;
         }
         double tpm = isoforms[i]._FPKM/sum_fpkm * 1e6;
         isoforms[i]._TPM = tpm;
         isoforms[i]._TPM_s = to_string(tpm);
      }

      for(auto it = isoforms.begin() ; it != isoforms.end();){
         if(it->_TPM < _kMinTPM || it->_FPKM < _kMinFPKM)
            it = isoforms.erase(it);
         else
            ++it;
      }
   }
//   else{
//      for(uint i=0; i< niso; ++i){
//         isoforms[i]._TPM = "FAILED";
//         isoforms[i]._FPKM = "FAILED";
//      }
//   }
   return success;
}

bool EmSolver::init(const int num_iso,
                  const vector<int> &count,
                  const vector<vector<double>> &model)
{
   int nrow = count.size();
   int ncol = num_iso;
   double total_count = accumulate(count.begin(), count.end(), 0.0);
   _theta = vector<double>(num_iso, total_count/num_iso);
   vector<size_t> erase_pos;
   for(int i= 0; i< nrow; ++i){
      bool remove = true;
      for(int j=0; j < ncol; ++j){
         if(model[i][j] > 1e-5) remove = false;
      }
      if(remove)
         erase_pos.push_back(i);
   }

   for(int i=0; i<nrow; ++i){
      if(binary_search(erase_pos.begin(), erase_pos.end(), i)) continue;
      _u.push_back(count[i]);
      _F.push_back(model[i]);
   }
   if(_u.empty()) return false;
//#ifdef DEBUG
//   nrow = _u.size();
//   Eigen::VectorXi obs_d(nrow);
//   Eigen::MatrixXd F(nrow, ncol);
//    for(size_t i =0; i < nrow; ++i){
//         for(size_t j= 0; j < ncol; ++j){
//            F(i,j) = _F[i][j];
//         }
//    }
//    for(size_t i = 0; i < nrow; ++i)
//      obs_d(i) = _u[i];
//    cerr<<"------n_i----------"<<endl;
//    cerr<<obs_d<<endl;
//    cerr<<"------Bias---------"<<endl;
//    cerr<<F<<endl;
//    cerr<<"-------------------"<<endl;
//#endif
    return true;
}

bool EmSolver::init( const int num_iso,
         const vector<int> &count,
         const vector<vector<double>> &model,
         const vector<vector<double>> &bias)
{
   int nrow = count.size();
   int ncol = num_iso;
   double total_count = accumulate(count.begin(), count.end(), 0.0);
   _theta = vector<double>(num_iso, total_count/num_iso);
   vector<size_t> erase_pos;
   for(int i= 0; i< nrow; ++i){
      bool remove = true;
      for(int j=0; j < ncol; ++j){
         if(model[i][j] > 1e-5) remove = false;
      }
      if(remove)
         erase_pos.push_back(i);
   }

   for(int i=0; i<nrow; ++i){
      if(binary_search(erase_pos.begin(), erase_pos.end(), i)) continue;
      _u.push_back(count[i]);
      _F.push_back(model[i]);
      _B.push_back(bias[i]);
   }
   if(_u.empty()) return false;

#ifdef DEBUG
   nrow = _u.size();
   Eigen::MatrixXd B(_B.size(), _B[0].size());
   Eigen::VectorXi obs_d(nrow);
   Eigen::MatrixXd F(nrow, ncol);

    for(size_t i =0; i < nrow; ++i){
         for(size_t j= 0; j < ncol; ++j){
            F(i,j) = _F[i][j];
         }
    }
    for(size_t i = 0; i < nrow; ++i)
      obs_d(i) = _u[i];
    for(size_t i =0; i< _B.size(); ++i){
       for(size_t j = 0; j< _B[0].size(); ++j)
          B(i,j) = _B[i][j];
    }
    cerr<<"------n_i----------"<<endl;
    cerr<<obs_d<<endl;
    cerr<<"------Sampleing Rates---------"<<endl;
    cerr<<F<<endl;
    cerr<<"------Bias-------------"<<endl;
    cerr<<B<<endl;
    cerr<<"------------------------"<<endl;
#endif
    return true;
}


bool EmSolver::run(){
   if (_u.empty()) return false;
   size_t nrow = _u.size();
   size_t ncol = _theta.size();

   Eigen::VectorXi obs_d(nrow);
   for(size_t i = 0; i < nrow; ++i)
      obs_d(i) = _u[i];

   Eigen::MatrixXd F(nrow, ncol);
     F.setZero(nrow, ncol);
   for(size_t i =0; i < nrow; ++i){ // updating the F matrix since theta changes
      for(size_t j= 0; j < ncol; ++j){
         F(i,j) = _F[i][j];
      }
   }

   Eigen::VectorXd theta(ncol);
   theta.setZero(ncol);
   for(int i = 0; i<ncol; ++i){
      theta[i] = _theta[i];
   }
   Eigen::VectorXd next_theta(ncol);
   next_theta.setZero(ncol);

   if(!_B.empty() && !_B[0].empty()){
      Eigen::MatrixXd B(_B.size(), _B[0].size());
      Eigen::VectorXd bias(_B[0].size());
      for(size_t i =0; i< _B.size(); ++i){
         for(size_t j = 0; j< _B[0].size(); ++j)
            B(i,j) = _B[i][j];
      }
      bias = (B.transpose() * B).inverse() * B.transpose() * Eigen::MatrixXd::Ones(_B.size(),1);
      Eigen::VectorXd next_bias(B.cols());
      next_bias.setZero(B.cols());

      for(int outer_it_num = 0; outer_it_num != _max_out_it_num; ++outer_it_num){

         for(int bias_it_num = 0; bias_it_num != _max_bias_it_num; ++bias_it_num){
            Eigen::MatrixXd U(B.rows(), B.cols());
            for(uint i = 0; i != B.rows(); ++i){
               double denom = bias.dot(B.row(i));
               if(denom < TOLERANCE) return false;
               for(uint j = 0; j != B.cols(); ++j)
                  U(i,j) = obs_d[i] * B(i,j)*bias[j]/denom;
            }
            for(uint j = 0; j != B.cols(); ++j){
               double numer = U.col(j).sum();
               double denom = 0.0;
               for(uint i = 0; i != B.rows(); ++i){
                  denom += B(i,j)* theta.dot(F.row(i));
               }
               if(denom < TOLERANCE) return false;
               next_bias[j] = numer/denom;
            }
            int num_changes = 0;
            for(uint j = 0; j < B.cols(); ++j){
               if( fabs( next_bias[j]-bias[j] )/next_bias[j] > _bias_change_limit ){
                  ++num_changes;
               }
               bias[j] = next_bias[j];
               next_bias[j] = 0.0;
            }
            if(num_changes == 0) break;
         }

         for(int theta_it_num = 0; theta_it_num < _max_theta_it_num; ++ theta_it_num){
            Eigen::MatrixXd U(nrow, ncol);
            for(uint i=0; i<nrow; ++i){
               double denom = theta.dot(F.row(i));
               if(denom < TOLERANCE) return false;
               for(uint j=0; j< ncol; ++j)
                  U(i,j) = obs_d[i] * F(i,j) * theta[j] / denom;
            }
            for(uint j=0; j<ncol; ++j){
               double numer = F.col(j).sum();
               double denom = 0;
               for (uint i =0; i< nrow; ++i){
                  denom = denom + F(i,j)*bias.dot(B.row(i));
               }
               if(denom < TOLERANCE) return false;
               next_theta[j] = numer / denom;
            }
            int num_changes = 0;
            for(uint j= 0; j< ncol; ++j){
               if( fabs(next_theta[j] - theta[j] )/next_theta[j] > _theta_change_limit )
                  ++num_changes;
               theta[j] = next_theta[j];
               next_theta[j] = 0.0;
            }
            if(num_changes == 0) break;
         }

         int num_changes = 0;
         for(uint j = 0; j< ncol; ++j){
            if( fabs(next_theta[j] - theta[j] )/next_theta[j] > _theta_change_limit )
                  ++num_changes;
         }
         for(uint j = 0; j < B.cols(); ++j){
               if( fabs( next_bias[j]-bias[j] )/next_bias[j] > _bias_change_limit ){
                  ++num_changes;
               }
         }
         if(num_changes == 0) break;
      }

      return true;
   }// end if


   for(int it_num = 0; it_num < _max_iter_num; ++it_num){
      /*E-step*/

      Eigen::MatrixXd U(nrow, ncol);
      U.setZero(nrow, ncol);

      for(size_t i =0; i < nrow; ++i){ // updating the unobserved data matrix
         double denom = F.row(i).dot(theta);
         if(denom < TOLERANCE) return false;
         for(size_t j= 0; j < ncol; ++j){
            double num = obs_d(i) * F(i,j)*theta[j];
            U(i,j) = num/denom;
         }
      }

      /*M-step*/
      for(int j = 0; j< U.cols(); ++j){
         //double den = F.col(j).sum();
         //if(den < TOLERANCE) return false;
         next_theta[j] = U.col(j).sum();
         //normalized_count += next_theta[j];
      }

      int num_changed = 0;
      for(size_t j=0; j<ncol; ++j){
         if((fabs(next_theta[j]-theta[j])/next_theta[j]) > _theta_change_limit)
            num_changed++;
         theta[j] = next_theta[j];
         next_theta[j] = 0.0;
      }

      if(num_changed == 0) {
         for(size_t j = 0; j< ncol; ++j){
         _theta[j] = theta[j];
         //cerr<<"isoform "<<j+1<<"'s raw read count: "<<_theta[j]<<endl;
         }
         break;
      }
   }

   return true;
}






