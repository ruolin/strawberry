/*
 * qp.cpp
 *
 *  Created on: Oct 15, 2015
 *      Author: ruolin
 */



#include "include/estimate.h"

#include <iostream>
#include <cassert>
#include <iterator>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <Eigen/Dense>
#include <stdexcept>
#include "contig.h"
// choose exact integral type

// program and solution types



uint generate_pair_end(const Contig& ct, const uint& start, int span,  SingleOrit_t orit){
/*
 * Given a singleton pair, its position 1bp before gap and orientation
 * infer the farthermost position of the other end.
 */
   vector<GenomicFeature> exons;
   for(auto gf:ct._genomic_feats){
      if(gf._match_op._code == Match_t::S_MATCH){
         exons.push_back(gf);
      }
   }
   if(orit == SingleOrit_t::Forward){
      auto it = lower_bound(exons.cbegin(), exons.cend(), start,
         [](const GenomicFeature &gf, uint s){return gf.right() < s;});
      //cout<<start<<" c "<< exons.cbegin()->left()<<endl;
      assert(it != exons.cend());
      span = span + start  - it->right();
      while(span >=0 && ++it != exons.cend()){
         span = span + it->left() - it->right();
      }
      if(it == exons.cend())
         return exons.crbegin()->right();
      return span + it->right();
   }


   if(orit == SingleOrit_t::Reverse){
      auto it = lower_bound(exons.crbegin(), exons.crend(), start,
            [](const GenomicFeature &gf, uint s){return gf.left() > s;});
      assert(it != exons.crend());
      span = span - (start - it->left());
      while(span >= 0 && ++it != exons.crend()){
         span = span - (it->right() - it->left());
      }
      if(it == exons.crend())
         return exons.cbegin()->left();
      return it->left() - span;
   }
   return 0;
}


ExonBin::ExonBin(std::set<uint> coordinate): _coords(coordinate){}

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

uint ExonBin::left_most() const
{
   return *_coords.begin();
}


uint ExonBin::right_most() const
{
   return *_coords.rbegin();
}


void ExonBin::add_read_mass(float mass)
{
   _whole_read_mass += mass;
}

bool ExonBin::add_frag(const Contig* fg)
{
   _frags.insert(fg);
   return true;
}

float ExonBin::read_count() const
{
   return _frags.size();
}

int ExonBin::bin_len() const
{
   int bin_len = 0;
   for(auto c = _coords.cbegin(); c != _coords.cend(); ++c){
      uint h = *(c++);
      bin_len += *c - h +1;
   }
   return bin_len;
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


Isoform::Isoform(Contig contig, int gene, int iso_id):
      _contig(contig), _gene_id(gene), _isoform_id(iso_id)
{
   _bais_factor = 0.0;
}

Estimation::Estimation(shared_ptr<InsertSize> insert_size, int read_len):
   _insert_size_dist(insert_size), _read_len(read_len) {}

void Estimation::set_maps(
             const int& iso_id,
             const int& fg_len,
             const float& mass,
             const Contig* read,
             const set<uint>& coords,
             map<set<uint>, ExonBin> & exon_bin_map,
             map<int, set<set<uint>>> &iso_2_bins_map)
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
      set<set<uint>> exon_bins = {coords};
      iso_2_bins_map.emplace(iso_id, exon_bins);
   }
   else{
      find_it2->second.insert(coords);
   }
}


void Estimation::overlap_exons(const vector<GenomicFeature> exons,
                     const uint left,
                     const uint right,
                     set<uint> &coords)
{
   for(auto const& gfeat : exons){
         if (GenomicFeature::overlap_in_genome(gfeat, left, right)){
            auto ret_l = coords.insert(gfeat.left());
            auto ret_r = coords.insert(gfeat.right());
            assert(ret_l.second && ret_r.second);
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
      map<set<uint>, ExonBin> & exon_bin_map,
      map<int, set<set<uint>>> &iso_2_bins_map)
{
   for(auto mp = hits.cbegin(); mp != hits.cend(); ++mp){

      map<set<uint>, pair<set<int>, int>> frag_mult_exonbin;
      double sr_fg_len = 0.0;
      if(mp->is_single_read() && infer_the_other_end){
         boost::mt19937 rng;
         double mean = _insert_size_dist->_mean;
         boost::normal_distribution<> nd(mean, _insert_size_dist->_sd);
         boost::variate_generator<boost::mt19937&, boost::normal_distribution<>> normal_rng(rng,nd);

         while( (sr_fg_len = normal_rng()) <= 0){}
      }

      for(auto iso = transcripts.cbegin(); iso != transcripts.cend(); ++iso){
         if(Contig::is_compatible(*mp, iso->_contig)){
            set<uint> coords;
            int frag_len = 0;

            /*For singleton, we random generate the other end */
            if(mp->is_single_read()){
               if(infer_the_other_end){
                  if(mp->single_read_orit() == SingleOrit_t::Forward){
                     uint other_end = generate_pair_end(iso->_contig, mp->right(), (int)sr_fg_len - _read_len, SingleOrit_t::Forward);
      //#ifdef DEBUG
      //                cout<<"Forward: "<<mp->left()<<":"<<other_end<<endl;
      //#endif
                     overlap_exons(exon_segs, mp->left(), other_end, coords);
                  }
                  else if (mp->single_read_orit() == SingleOrit_t::Reverse) {
                     uint other_end = generate_pair_end(iso->_contig, mp->left(), (int)sr_fg_len - _read_len, SingleOrit_t::Reverse);
      //#ifdef DEBUG
      //                  cout<<"Reverse: "<<other_end<<":"<<mp->right()<<endl;
      //#endif
                     overlap_exons(exon_segs, other_end, mp->right(), coords);
                  }
                  auto ret = frag_mult_exonbin.find(coords);
                  if(ret == frag_mult_exonbin.end()){
                     frag_mult_exonbin.emplace (coords, pair<set<int>,int>({iso->_isoform_id}, sr_fg_len));
                  }
                  else{
                     ret->second.first.insert(iso->_isoform_id);
                  }
               } // end if infer the other
               else{
                  overlap_exons(exon_segs, mp->left(), mp->right(), coords);
                  frag_len = Contig::exonic_overlaps_len(iso->_contig, mp->left(), mp->right());
                  set_maps(iso->_isoform_id, frag_len, mp->mass(), &(*mp), coords, exon_bin_map , iso_2_bins_map);
               }
            } // and and single end

            else{
               set<uint> l;
               overlap_exons(exon_segs,mp->left(), mp->gap_left()-1, l);
               for(auto i: l)
                  coords.insert(i);
               set<uint> r;
               overlap_exons(exon_segs, mp->gap_right()+1, mp->right(), r);
               for(auto i: r)
                  coords.insert(i);
               frag_len = Contig::exonic_overlaps_len(iso->_contig, mp->left(), mp->right());
               set_maps(iso->_isoform_id, frag_len, mp->mass(), &(*mp), coords, exon_bin_map , iso_2_bins_map);
            }
            //assert(!coords.empty());
//#ifdef DEBUG
//            cout<<"coords";
//            for(auto const& c : coords)
//               cout<<" "<<c;
//            cout<<" hit: "<<mp->left()<<endl;
//#endif
//            cout<<frag_len<<""<<endl;
         } // end if compatible condition
      }// end second inner for loop
      if(frag_mult_exonbin.size() == 1){
         const set<uint> coords  = frag_mult_exonbin.begin()->first;
         const int frag_len = frag_mult_exonbin.begin()->second.second;
         for(auto iso_num: frag_mult_exonbin.begin()->second.first)
            set_maps(iso_num, frag_len, mp->mass(), &(*mp),coords, exon_bin_map , iso_2_bins_map);
      }
      if(frag_mult_exonbin.size() > 1){
         for(auto iso = transcripts.cbegin(); iso != transcripts.cend(); ++iso){
            if(Contig::is_compatible(*mp, iso->_contig)){
               set<uint> coords;
               overlap_exons(exon_segs, mp->left(), mp->right(), coords);
               int frag_len = Contig::exonic_overlaps_len(iso->_contig, mp->left(), mp->right());
               set_maps(iso->_isoform_id, frag_len, mp->mass(), &(*mp), coords, exon_bin_map , iso_2_bins_map);
            }
         }
      }
   }
}

void Estimation::calculate_raw_iso_counts(const map<int, set<set<uint>>> &iso_2_bins_map,
          const map<set<uint>, ExonBin> & exon_bin_map)
{

}


void Estimation::calcuate_bin_weight( const map<int, set<set<uint>>> &iso_2_bins_map,
                                  const map<int,int> &iso_2_len_map,
                                  const int m,
                                  map<set<uint>, ExonBin> & exon_bin_map)
{
   //map<int, double> iso_2_bias_factor_map;
   for(auto it = iso_2_bins_map.cbegin(); it != iso_2_bins_map.cend(); ++it){
      for(auto const &bin_coord : it->second){
         double weight = 0.0;
         //iso_num = exon_bin_map.at(bin_coord)
         for(auto const &len_mass_pair: exon_bin_map.at(bin_coord)._iso_2_frag_lens.at(it->first)){
            int num_start_pos = iso_2_len_map.at(it->first)-len_mass_pair.first+1;
            weight += len_mass_pair.second * _insert_size_dist->emp_dist_pdf(len_mass_pair.first) / num_start_pos;
         }
         exon_bin_map.at(bin_coord)._iso_bias_map[it->first] = m*weight;
      }
      //iso_2_bias_factor_map[it->first] = bias_iso;
      //cout<<"bias iso "<<bias_iso<<endl;
   }

   for(auto & il: iso_2_len_map){
#ifdef DEBUG
            cout<<"isoform length for "<<il.first<<": "<<il.second<<endl;
#endif
   }

//#ifdef DEBUG
//   for(auto const & bin: exon_bin_map){
//      cout<<"bin_coord: ";
//      for(auto it = bin.first.cbegin() ; it!= bin.first.cend(); ++it)
//         cout<<*it<<"-";
//      cout<<endl;
//      for(auto f : bin.second._frags)
//         cout<<"read "<<f->left()<<"-"<<f->right()<<"\t";
//      cout<<endl;
////      for(auto const & iso: bin.second._iso_bias_map){
////         cout<<"iso bias: "<<iso.first<<": "<<iso.second<<endl;;
////      }
//   }
//#endif
}

bool Estimation::estimate_abundances(const map<set<uint>, ExonBin> & exon_bin_map,
                     const double mass,
                     vector<Isoform>& isoforms,
                     bool use_qp,
                     bool with_bias_correction)
{
   size_t nrow = exon_bin_map.size();
   size_t niso = isoforms.size();

   vector<int> n (nrow);
   vector<vector<double>> alpha (nrow, vector<double>(niso));
   unsigned i = 0;
   for(auto const& bin: exon_bin_map){
      n[i] = bin.second.read_count();
      for(unsigned j = 0; j< niso; ++j){
         auto ret = bin.second._iso_bias_map.find(j+1);
         if (ret == bin.second._iso_bias_map.end())
            alpha[i][j] = 0;
         else
            alpha[i][j] = ret->second;
      } // //inner loop
      ++i;
   } // outer loop
   EmSolver em(niso, n, alpha);
   bool success = em.run();
   if(success){
      for(int i=0; i< niso; ++i)
         isoforms[i]._abundance = to_string(em._theta[i]);
   }
   else{
      for(int i=0; i< niso; ++i)
         isoforms[i]._abundance = "FAILED";
   }
   return success;

      /*
       * Now using EM algorithm
       */
//   if(use_qp){
//      /*
//       * Set objective function for CGAL. x'Dx + c'x + c_0
//       * Quadratic term -> linear term -> constant term
//       */
//      Program qp (CGAL::EQUAL, true, 0, false, 0);
//      matrix D = 2*boost::numeric::ublas::prod(trans(alpha), alpha);
//      double c_0 = boost::numeric::ublas::inner_prod(b,b);
//      ublas_vector c = -2 *boost::numeric::ublas::prod(trans(alpha), b);
//      /*
//       * Set equaltiy constriant  1'x = 1;
//       */
//   //   for(unsigned idx = 0;  idx < ncol; ++idx){
//   //      qp.set_a(idx, 0, 1);
//   //   }
//   //   qp.set_b(0,1);
//
//      i = 0;
//      unsigned j = 0;
//      for(matrix::iterator1 it1 = D.begin1(); it1 != D.end1(); ++it1){
//         for(matrix::iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2){
//            if(j>i) continue;
//            qp.set_d(i, j, *it2);
//            ++j;
//         }
//         j = 0;
//         ++i;
//      }
//
//      i = 0;
//      for(ublas_vector::iterator it = c.begin(); it != c.end(); ++it){
//         qp.set_c(i, *it);
//         ++i;
//      }
//      qp.set_c0(c_0);
//
//      /*
//       * run solver
//       */
//      Solution s = CGAL::solve_nonnegative_quadratic_program(qp, ET());
//      assert (s.solves_nonnegative_quadratic_program(qp));
//
//   #ifdef DEBUG
//   //      if(s.objective_value() < 0){
//   //          CGAL::print_nonnegative_quadratic_program
//   //          (std::cout, qp, "first_nonnegative_qp");
//   //          cout<<s;
//   //      }
//   #endif
//      // output solution
//      //cout<<s;
//         i = 0;
//      for(auto it = s.variable_values_begin(); it != s.variable_values_end(); ++it){
//         isoforms[i]._abundance = it->num.to_double()/it->den.to_double();
//         //std::stringstream ss;
//         //ss<<*it;
//         cout<<isoforms[i]._abundance<<endl;
//         ++i;
//      }
//   }
//
//   else{
//   }
}

EmSolver::EmSolver(const int num_iso, const vector<int> &count, const vector<vector<double>> &model):
      _num_isoforms(num_iso),
      _theta(num_iso, 1.0/num_iso)
      //_u(count),
      //_F(model)
{
   int nrow = count.size();
   int ncol = _num_isoforms;
   vector<size_t> erase_pos;
   for(int i= 0; i< nrow; ++i){
      bool remove = true;
      for(int j=0; j < ncol; ++j){
         if(model[i][j] > 1e-6) remove = false;
      }
      if(remove)
         erase_pos.push_back(i);
   }

   for(int i=0; i<nrow; ++i){
      if(binary_search(erase_pos.begin(), erase_pos.end(), i)) continue;
      _u.push_back(count[i]);
      _F.push_back(model[i]);
   }
#ifdef DEBUG
   nrow = _u.size();
   Eigen::VectorXi obs_d(nrow);
   //ublas_vector obs_d (nrow);

   Eigen::MatrixXd F(nrow, ncol);
   //matrix F(nrow, ncol);
    for(size_t i =0; i < nrow; ++i){ // updating the F matrix since theta changes
         for(size_t j= 0; j < ncol; ++j){
            F(i,j) = _F[i][j];
         }
    }
    for(size_t i = 0; i < nrow; ++i)
      obs_d(i) = _u[i];
    cout<<"------n_i----------"<<endl;
    cout<<obs_d<<endl;
    cout<<"------Bias---------"<<endl;
    cout<<F<<endl;
    cout<<"-------------------"<<endl;
#endif
}


bool EmSolver::run(){
   if (_u.empty()) return false;

   size_t nrow = _u.size();
   size_t ncol = _num_isoforms;

   //ublas_vector obs_d (nrow);
   Eigen::VectorXi obs_d(nrow);
   for(size_t i = 0; i < nrow; ++i)
      obs_d(i) = _u[i];

   //matrix F_init(nrow, ncol);
   Eigen::MatrixXd F_init(nrow,ncol);
   for(size_t i =0; i < nrow; ++i)
         for(size_t j= 0; j < ncol; ++j)
            F_init(i,j) = _F[i][j];

   //ublas_vector next_theta(ncol,0.0);
   Eigen::VectorXd next_theta(ncol);
   next_theta.setZero(ncol);

   for(int it_num = 0; it_num < _max_iter_num; ++it_num){
      /*
       * E-step
       * */
      //matrix F(nrow, ncol, 0.0);
      //ublas_vector F_row_sum(nrow, 0.0);
      //matrix U (nrow, ncol, 0.0);
      Eigen::MatrixXd F(nrow, ncol);
      F.setZero(nrow, ncol);
      Eigen::VectorXd F_row_sum(nrow);
      F_row_sum.setZero(nrow);
      Eigen::MatrixXd U(nrow, ncol);
      U.setZero(nrow, ncol);

      for(size_t i =0; i < nrow; ++i){ // updating the F matrix since theta changes
         for(size_t j= 0; j < ncol; ++j){
            F(i,j) = _F[i][j] * _theta[j];
            F_row_sum(i) += F(i,j);
         }
      }
      for(size_t i =0; i < nrow; ++i){ // updating the unobserved data matrix
         for(size_t j= 0; j < ncol; ++j){
            double num = obs_d(i) * F(i,j);
            if(F_row_sum(i) < TOLERANCE) return false;
            double denom = F_row_sum(i);
            U(i,j) = num/denom;
         }
      }
      /*
       * M-step
       */
      double normalized_count = 0.0;
      for(size_t j = 0; j< U.cols(); ++j){
         //matrix_col U_col(U,j);
         //one_vector one(nrow, 1);
         //matrix_col F_init_col(F_init,j);

         //double denom = boost::numeric::ublas::inner_prod(one, F_init_col);
         //next_theta[j] = boost::numeric::ublas::inner_prod(one, U_col)/denom;
         next_theta[j] = U.col(j).sum();
         normalized_count += next_theta[j];
      }

      /* Normalizing next_theta */
      for(size_t j = 0; j< ncol; ++j){
         next_theta[j] /= normalized_count;
      }

      int num_changed = 0;
      for(size_t j=0; j<ncol; ++j){
         if(next_theta[j] > _theta_limit && (fabs(next_theta[j]-_theta[j])/next_theta[j]) > _theta_change_limit)
            num_changed++;
         _theta[j] = next_theta[j];
         next_theta[j] = 0.0;
      }

      if(num_changed == 0) {
         for(size_t j = 0; j< ncol; ++j){
         cout<<"tehta: "<<_theta[j]<<endl;
         }
         break;
      }
   }

   return true;
}







