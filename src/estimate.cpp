
#include <iostream>
#include <cassert>
#include <iterator>
#include <random>
#include <Eigen/Dense>
#include <stdexcept>
#include "estimate.hpp"
#include "fasta.h"
#include "contig.h"
#include "bias.hpp"
// choose exact integral type

// program and solution types
using namespace std;
const double LocusContext::_kMinFrac = kMinIsoformFrac;
//const double LocusContext::_kMinTPM = 0;

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


set<pair<uint,uint>> LocusContext::overlap_exons(const vector<GenomicFeature>& exons, const Contig& read) const
{
   set<pair<uint,uint>> coords;
   for(auto const& gfeat : exons){
      if (gfeat._match_op._code != Match_t::S_MATCH) continue;

      for(auto const &read_f: read._genomic_feats){
         if (read_f._match_op._code != Match_t::S_MATCH) continue;

         if (GenomicFeature::overlaps(read_f, gfeat)){
            auto ret = coords.insert(pair<uint,uint>(gfeat.left(),gfeat.right()));
            //assert(ret.second);
         }
      }
   }
   return coords;
}



void LocusContext::assign_exon_bin(
      const vector<Contig> &hits,
      const vector<GenomicFeature> & exon_segs)
{
/*
 * assign reads and transcripts to exon bin.
 */
   for(auto mp = hits.cbegin(); mp != hits.cend(); ++mp){

      double sr_fg_len = 0.0;
      if(mp->is_single_read() && infer_the_other_end){ // currently infer_the_other_end is always disabled
         //random_device rd;
         //mt19937 gen(rd());
         mt19937 gen(3); // we use a fixed seed to make sure the output are the same every time.
         double mean = _sample._insert_size_dist->_mean;
         double sd = _sample._insert_size_dist->_sd;
         normal_distribution<> nd(mean, sd);

         while( (sr_fg_len = nd(gen)) <= 0){}
      }

      for(auto iso = _transcripts.cbegin(); iso != _transcripts.cend(); ++iso){
         if(Contig::is_compatible(*mp, iso->_contig)){
            set<pair<uint,uint>> coords;
            int frag_len = 0;
            //Bias::iso_bias(*mp, *iso);
            /*For singleton, we random generate the other end */
            if(mp->is_single_read()){
//               if(infer_the_other_end){ // currently infer_the_other_end is always disabled
//                  Contig new_read;
//                  if(mp->single_read_orit() == SingleOrit_t::Forward){
//                     uint other_end = generate_pair_end(iso->_contig, *mp, (int)sr_fg_len - _read_len, SingleOrit_t::Forward, new_read);
//                     if(other_end == 0) continue;
//                     coords = overlap_exons(exon_segs, new_read);
//                  }
//
//                  else if (mp->single_read_orit() == SingleOrit_t::Reverse) {
//
//                     uint other_end = generate_pair_end(iso->_contig, *mp, (int)sr_fg_len - _read_len, SingleOrit_t::Reverse, new_read);
//                     if(other_end == 0) continue;
//                     assert(coords.empty());
//                     coords = overlap_exons(exon_segs, new_read);
//                  }
//                  frag_len = Contig::exonic_overlaps_len(iso->_contig, new_read.left(), new_read.right());
//                  set_maps(iso->id(), frag_len, new_read.mass(), new_read, coords);
//               } // end if infer the other

               //else{
                  coords = overlap_exons(exon_segs, *mp);
                  frag_len = Contig::exonic_overlaps_len(iso->_contig, mp->left(), mp->right());
                  set_maps(iso->id(), frag_len, mp->mass(), *mp, coords);
               //}
            } // and and single end

            else{
               coords = overlap_exons(exon_segs,*mp);
               frag_len = Contig::exonic_overlaps_len(iso->_contig, mp->left(), mp->right());
               set_maps(iso->id(), frag_len, mp->mass(), *mp, coords);
            }

         } // end if compatible condition
      }// end second inner for loop
   }
}


void LocusContext::set_theory_bin_weight() {

   for(auto it = iso_2_bins_map.cbegin(); it != iso_2_bins_map.cend(); ++it){
      assert(_transcripts[it->first].id() == it->first);
      for(auto const &bin_idx: it->second){
         double weight = 0.0;
         vector<pair<uint,uint>> exon_segs;
         vector<uint> implicit_exon_idx = exon_bins.at(bin_idx).bin_under_iso(_transcripts[it->first], exon_segs);
         vector<uint> seg_lens(exon_segs.size());
         for(uint i=0; i< seg_lens.size(); ++i){
            seg_lens[i] = exon_segs[i].second - exon_segs[i].first + 1;
         }

         int lmax = accumulate(seg_lens.begin(), seg_lens.end(),0);
         int lmin;
         if(_sample._insert_size_dist->_use_emp)
            lmin = _sample._insert_size_dist->_start_offset;
         else
            //lmin = 2*_read_len;
            lmin = _read_len;
         if(seg_lens.size() > 2)
            lmin = max(lmin, accumulate(seg_lens.begin() +1, seg_lens.end()-1, 0));
         //std::cerr<<"min "<<lmin<<" max "<<lmax<<std::endl;
         for(int fl = lmin ; fl <= lmax; ++fl){
            double le_eff = exon_bins.at(bin_idx).effective_len(seg_lens,implicit_exon_idx, fl, _read_len);
            double tmp = _sample._insert_size_dist->emp_dist_pdf(fl)* le_eff / (_transcripts[it->first]._length - fl + 1);
            weight += tmp;
         }
         //std::cerr<<"weight "<<weight<<std::endl;
         exon_bins.at(bin_idx)._bin_weight_map[it->first] = weight;

      }
   }
}

void LocusContext::set_bin_weight_without_frag_dist() {
   for(auto it = iso_2_bins_map.cbegin(); it != iso_2_bins_map.cend(); ++it){
      assert(_transcripts[it->first].id() == it->first);
      for(auto const &bin_idx: it->second){
         double weight = 1.0 / _transcripts[it->first]._length;
         //double weight = 1.0 ;
         exon_bins.at(bin_idx)._bin_weight_map[it->first] = weight;

      }
   }

}


void LocusContext::set_empirical_bin_weight(const map<int,int> &iso_2_len_map, const int m) {
   map<int, double> iso_weight_factor_map;
   for(auto it = iso_2_bins_map.cbegin(); it != iso_2_bins_map.cend(); ++it){
      double iso_weight = 0;
      for(auto const &bin_idx : it->second){
         double weight = 0.0;
         for(auto const &len_mass_pair: exon_bins.at(bin_idx)._iso_2_frag_lens.at(it->first)){
            int num_start_pos = iso_2_len_map.at(it->first)-len_mass_pair.first+1;
            weight += len_mass_pair.second * _sample._insert_size_dist->emp_dist_pdf(len_mass_pair.first) / num_start_pos;
         }
         exon_bins.at(bin_idx)._bin_weight_map[it->first] = weight;
         iso_weight += weight;
      }
      iso_weight_factor_map[it->first] = iso_weight;
   }

   for(auto it = exon_bins.begin(); it != exon_bins.end(); ++it){
      for(auto iso = iso_weight_factor_map.cbegin(); iso != iso_weight_factor_map.cend(); ++iso){
         it->_bin_weight_map[iso->first]  /= iso->second;
      }
   }

   for(auto & il: iso_2_len_map){
//#ifdef DEBUG
//            cout<<"isoform length for "<<il.first<<": "<<il.second<<endl;
//#endif
   }
}

bool LocusContext::estimate_abundances()
{
   size_t nrow = exon_bins.size();
   size_t niso = _transcripts.size();

   vector<int> n (nrow);
   vector<vector<double>> alpha(nrow, vector<double>(niso));
   unsigned i = 0;
   for(auto const& bin: exon_bins){
      n[i] = bin.read_count();
      for(unsigned j = 0; j< niso; ++j){
         auto ret = bin._bin_weight_map.find(j);
         if (ret == bin._bin_weight_map.end())
            alpha[i][j] = 0.0;
         else
            alpha[i][j] = ret->second;
      } // //inner loop
      ++i;
   } // outer loop
#ifdef DEBUG
   for(uint i=0; i!= alpha.size();++i){
      for(uint j=0; j!= alpha[i].size(); ++j)
         cout<<alpha[i][j]<<" ";
      cout<<endl;
   }
#endif
   bool success;
   EmSolver em;
   success = em.init(niso, n, alpha);
   if(success) em.run();

   if(success){
      for(uint i=0; i<niso; ++i){
         fprintf(_p_log_file, "isoform %d has %f raw read count.\n", i+1, em._theta[i]);
      }
      double sum_fpkm = 0.0;
      for(uint i=0; i< niso; ++i){
         double kb = 0.0;
         if(effective_len_norm){
            kb = _transcripts[i]._length - _sample._insert_size_dist->_mean;
            if (kb < 0){
               _transcripts[i]._FPKM_s = "nan";
               continue;
            }
            kb =1e3/kb;
         }
         else{
            kb = 1e3/_transcripts[i]._length;
         }
         double rpm = 1e6/_sample.total_mapped_reads();
         double fpkm = em._theta[i]*rpm*kb;
//#ifdef DEBUG
//         cout<<"theta: "<<i<<" exp: "<<em._theta[i]<<endl;
//#endif
         _transcripts[i]._FPKM = fpkm;
         sum_fpkm += fpkm;
         _transcripts[i]._FPKM_s = to_string(fpkm);
      }
      for(uint i=0; i< niso; ++i){
         if(_transcripts[i]._FPKM_s == "nan"){
            _transcripts[i]._frac_s = "nan";
            continue;
         }
         double frac = _transcripts[i]._FPKM/sum_fpkm ;
         _transcripts[i]._frac = frac;
         _transcripts[i]._frac_s = to_string(frac);
      }
      if(filter_by_expression){
         for(auto it = _transcripts.begin() ; it != _transcripts.end();){
            if(it->_frac < kMinIsoformFrac) {
               it = _transcripts.erase(it);
            }
            else {
               ++it;
            }
         }
      }
   }
//   else{
//      for(uint i=0; i< niso; ++i){
//         _transcripts[i]._TPM = "FAILED";
//         _transcripts[i]._FPKM = "FAILED";
//      }
//   }
   return success;
}

bool EmSolver::init(const int num_iso,
                  const vector<int> &count,
                  const vector<vector<double>> &model)
{
   //std::cerr << "count size: " << count;
   //std::cerr << "model: " << model << std::endl;
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
#ifdef DEBUG
   nrow = _u.size();
   Eigen::VectorXi obs_d(nrow);
   Eigen::MatrixXd F(nrow, ncol);
    for(size_t i =0; i < nrow; ++i){
         for(size_t j= 0; j < ncol; ++j){
            F(i,j) = _F[i][j];
         }
    }
    for(size_t i = 0; i < nrow; ++i)
      obs_d(i) = _u[i];
    cerr<<"------n_i----------"<<endl;
    cerr<<obs_d<<endl;
    cerr<<"-------------------"<<endl;
    cerr<<F<<endl;
#endif
    return true;
}

bool EmSolver::run(){
   if (_u.empty()) return false;
   size_t nrow = _u.size();
   size_t ncol = _theta.size();

   //initialization
   Eigen::VectorXi obs_d(nrow);
   for(size_t i = 0; i < nrow; ++i)
      obs_d(i) = _u[i];

   // Conditional probability update matrix
   Eigen::MatrixXd F(nrow, ncol);
   F.setZero(nrow, ncol);
   for(size_t i =0; i < nrow; ++i){
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

   Eigen::MatrixXd U(nrow, ncol);
   U.setZero(nrow, ncol);

   Eigen::MatrixXd newF(nrow, ncol);
   newF.setZero(nrow, ncol);

   for(int it_num = 0; it_num < _max_iter_num; ++it_num){ // if no bias est
      /*E-step*/

      // posterior latent class probability

      for(size_t i =0; i < nrow; ++i){ // updating the unobserved data matrix
         double denom = F.row(i).dot(theta);
         if(denom == 0) {
            return false;
         }
         for(size_t j= 0; j < ncol; ++j){
            double num = obs_d(i) * F(i,j)*theta[j];
            U(i,j) = num/denom;
         }
      }


      /*M-step*/
      for(int j = 0; j< U.cols(); ++j){
         next_theta[j] = U.col(j).sum();
      }

      for(size_t j= 0; j < ncol; ++j){
         double denom = F.col(j).sum();
         for(size_t i =0; i < nrow; ++i){
            if(denom == 0) {
               newF(i,j) == 0;
               //return false;
            } else {
               newF(i, j) = F(i, j) / denom;
            }
         }
      }

      F = newF;
      Eigen::VectorXd dist = next_theta - theta;
      if( dist.norm()  < _theta_change_limit) break;
      theta = next_theta;
   }

   for(size_t j = 0; j< ncol; ++j){
      _theta[j] = theta[j];
   }
   return true;
}

