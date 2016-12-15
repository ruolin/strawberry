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

#include <algorithm>
#include <iterator>
#include <random>
#include <math.h>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include "alignments.h"
#include "fasta.h"
#include "assembly.h"
#include "estimate.hpp"
#include "bias.hpp"
#include "interval.hpp"
#include "utils.h"

using namespace std;


#if ENABLE_THREADS
mutex out_file_lock;
mutex thread_pool_lock;
atomic<int> curr_thread_num = {0};

void decr_pool_count()
{
   --curr_thread_num;
}
#endif
/*
 * Global utility function begin:
 */
bool hit_lt_cluster(const ReadHit& hit, const HitCluster& cluster, uint olap_radius){
   if(hit.ref_id() != cluster.ref_id())
     return hit.ref_id() < cluster.ref_id();
   else
     return hit.right() + olap_radius < cluster.left();
}

bool hit_gt_cluster(const ReadHit& hit, const HitCluster& cluster, uint olap_radius){

   if(hit.ref_id() != cluster.ref_id()){
     //cout<<"shouldn't\t"<<hit.ref_id()<<":"<<cluster.ref_id()<<endl;
     return hit.ref_id() > cluster.ref_id();
   }
   else{
     return hit.left() > cluster.right() + olap_radius;
   }
}

bool hit_complete_within_cluster(const PairedHit& hit, HitCluster& cluster, uint olap_radius){
   if(hit.ref_id() != cluster.ref_id()){
     return false;
   }
   else{
     if(hit.right_pos()+olap_radius < cluster.left())
       return false;
     if(cluster.right()+olap_radius < hit.left_pos())
       return false;
   }
   return true;
}

/*
 * Global utility functions end:
 */

HitCluster::HitCluster():
   _leftmost(UINT_MAX),
   _rightmost(0),
   _plus_strand_num_hits(0),
   _minus_strand_num_hits(0),
   _first_encounter_strand(Strand_t::StrandUnknown),
   _ref_id(-1),
   _final(false),
   _raw_mass(0.0)
{}


RefID HitCluster::ref_id() const
{
   return _ref_id;
}

void HitCluster::ref_id(RefID id)
{
   _ref_id = id;
}

void HitCluster::addRefContig(Contig *contig)
{
   if(_ref_id != -1){
     assert(_ref_id == contig->ref_id());
   }
   else{
     _ref_id = contig->ref_id();
   }
   if (gene_id().empty()) gene_id() = contig->parent_id();
   else assert(gene_id() == contig->parent_id());
   _leftmost = min(_leftmost, contig->left());
   _rightmost = max(_rightmost, contig->right());
   _ref_mRNAs.push_back(contig);
}

uint HitCluster::left() const
{
   return _leftmost;
}
uint HitCluster::right() const
{
   return _rightmost;
}

void HitCluster::left(uint left)
{
   _leftmost = left;
}
void HitCluster::right(uint right)
{
   _rightmost = right;
}

int HitCluster::size() const {
   return _hits.size() + _open_mates.size();
}

int HitCluster::len() const{
   return _rightmost - _leftmost + 1;
}

Strand_t HitCluster::ref_strand() const{
   assert(!_ref_mRNAs.empty());
   Strand_t strand = _ref_mRNAs[0]->strand();
   for(auto &i: _ref_mRNAs){
     assert(i->strand() == strand);
   }
   return strand;
}

Strand_t HitCluster::guessStrand() const{
   if(_first_encounter_strand == Strand_t::StrandUnknown){
     LOG("HitCluster ", _ref_id, ":", _leftmost,"-",_rightmost," does not have strand information. It is likely to be a single exon transcript." );
     return Strand_t::StrandUnknown;
   }
   if(_first_encounter_strand == Strand_t::StrandPlus){
     if(_plus_strand_num_hits >=3 ) return Strand_t::StrandPlus;
     else{
       if(_plus_strand_num_hits >= _minus_strand_num_hits) return Strand_t::StrandPlus;
       else return Strand_t::StrandMinus;
     }
   }
   else{
     if(_minus_strand_num_hits >=3) return Strand_t::StrandMinus;
     else{
       if(_minus_strand_num_hits >= _plus_strand_num_hits) return Strand_t::StrandMinus;
       else return Strand_t::StrandPlus;
     }
   }
}

void HitCluster::reset(uint old_left,
                  uint old_right,
                  RefID old_ref_id,
                  int old_plus_hit,
                  int old_minus_hit,
                  Strand_t old_strand,
                  double old_mass){
   _leftmost = old_left;
   _rightmost = old_right;
   _ref_id = old_ref_id;
   _plus_strand_num_hits = old_plus_hit;
   _minus_strand_num_hits = old_minus_hit;
   _first_encounter_strand = old_strand;
   _raw_mass = old_mass;
}

void HitCluster::reset(uint old_left,
                  uint old_right,
                  RefID old_ref_id
                  ){
   _leftmost = old_left;
   _rightmost = old_right;
   _ref_id = old_ref_id;
}


bool HitCluster::addHit(const PairedHit &hit){

   if(_final){
     return false;
   }
   assert(_ref_id == hit.ref_id());
   if(hit.contains_splice()){
     if(hit.strand() == Strand_t::StrandPlus){
       if(_minus_strand_num_hits == 0 && _plus_strand_num_hits == 0){
         _first_encounter_strand = Strand_t::StrandPlus;
       }
       ++_plus_strand_num_hits;
     }
     else if(hit.strand() == Strand_t::StrandMinus){
       if(_minus_strand_num_hits == 0 && _plus_strand_num_hits == 0){
         _first_encounter_strand = Strand_t::StrandMinus;
       }
       ++_minus_strand_num_hits;
     }
     else{
       //This cause problem with some versions of STAR
       //assert(false);
     }
   }
#ifdef DEBUG
   if(hit._left_read){
      assert(hit.left_read_obj().cigar().front()._type == MATCH ||
          hit.left_read_obj().cigar().front()._type == SOFT_CLIP);
   }
   if(hit._right_read){
       assert(hit.right_read_obj().cigar().back()._type == MATCH||
            hit.right_read_obj().cigar().back()._type == SOFT_CLIP);
   }
#endif
   _hits.push_back(hit);
   return true;
}

//void HitCluster::count_current_intron(const std::vector<std::pair<uint,uint>> &intron)
//{
//   bool not_found = true;
//   for(auto it: intron){
//     auto p = _current_intron_counter.find(it);
//     if(p != _current_intron_counter.end()) {
//       not_found = false;
//       p->second++;
//     }
//   }
//   if(not_found){
//     _current_intron_counter.clear();
//     for(auto it: intron){
//       _current_intron_counter.emplace(it,0);
//     }
//   }
//}
//
//bool HitCluster::current_intron_is_reliable() const
//{
//   for(auto it:_current_intron_counter){
//     if(it.second > kMinJuncSupport)
//       return true;
//   }
//   return false;
//}

bool HitCluster::addOpenHit(const ReadHitPtr hit, bool extend_by_hit, bool extend_by_partner)
{
   if(this->size() > kMaxFragPerCluster)
     return false;
   uint orig_left = _leftmost;
   uint orig_right = _rightmost;
   RefID orig_ref_id = _ref_id;
   uint hit_left = hit->left();
   uint hit_right = hit->right();
   Strand_t hit_strand = hit->strand();
   RefID hit_ref_id = hit->ref_id();
   uint hit_partner_pos = hit->partner_pos();
   ReadID hit_id = hit->read_id();

//   if(hit->contains_splice() && hit->intron_len() > kMaxIntronLen4ExtCluster){
//     vector<pair<uint,uint>> intron_coords = hit->intron_coords();
//
//     //count_current_intron(intron_coords);
//     if(extend_by_hit && current_intron_is_reliable()){
//       _leftmost = min(_leftmost, hit_left);
//       _rightmost = max(_rightmost, hit_right);
//     }
//     if(extend_by_partner && hit_partner_pos != 0 && current_intron_is_reliable()){
//       if((int)hit_partner_pos - (int)hit_left > max_inner_d){
//         LOG_ERR("Read Pair ", hit_left, "-",hit_partner_pos, " inner distance is larger than ",  max_inner_d);
//         return false;
//       }
//       _rightmost = max(max(_rightmost, hit->right()), hit->partner_pos());
//     }
//   } //end if
//   else{
     if(extend_by_hit){
       _leftmost = min(_leftmost, hit_left);
       _rightmost = max(_rightmost, hit_right);
     }
     if(extend_by_partner && hit_partner_pos != 0 && hit->partner_ref_id() == _ref_id){
       if((int)hit_partner_pos - (int)hit_left < kMaxIntronLength){
         _rightmost = max(max(_rightmost, hit->right()), hit->partner_pos());
       }
     }
//   }//end else


   // Double check. This is only useful when called in Sample::mergeCluster()
   if(hit_lt_cluster(*hit, *this, kMaxOlapDist)){
     reset(orig_left, orig_right, orig_ref_id);
     return false;
   }

   if(abs((int)hit_right - (int)hit_left) > kMaxFragSpan){
     reset(orig_left, orig_right, orig_ref_id);
     LOG_WARN("Hit start at ",hit_left, "  is longer than max gene length, skipping");
     return false;
   }

   if(_ref_id == -1){
     if(hit_ref_id != -1)
       _ref_id = hit_ref_id;
   } else{
     assert(_ref_id == hit_ref_id);
   }

   if(hit->is_singleton() || hit->partner_ref_id() != _ref_id){
     if(hit->reverseCompl()){
       PairedHit ph(nullptr, hit);
       addHit(ph);
     }
     else{
       PairedHit ph(hit, nullptr);
       addHit(ph);
     }
   }

   else{
     unordered_map<ReadID, list<PairedHit>>::iterator iter_open = _open_mates.find(hit_id);
     if( iter_open == _open_mates.end()){


       if(hit->partner_pos() > hit->right()){
         if(hit->reverseCompl())
         {
            LOG_ERR("Possible wrong read orientation at chr: ", hit->ref_id(), " for read start at ", hit->left(), " and his partner at ",hit->partner_pos() );
         }
         PairedHit open_hit(hit, nullptr);
         unordered_map<ReadID, list<PairedHit>>::iterator ins_pos;
         list<PairedHit> chain;
         chain.push_back(move(open_hit));
         bool status;
         tie(ins_pos, status) = _open_mates.insert(make_pair(hit_id, chain));
         assert(status);
       }
       else if(hit->partner_pos() < hit->left()){
         if(!hit->reverseCompl())
         {
            LOG_ERR("Possible wrong read orientation at chr: ", hit->ref_id(), " for read start at ", hit->left(), " and his partner at ",hit->partner_pos() );
         }

         PairedHit open_hit(nullptr, hit);
         unordered_map<ReadID, list<PairedHit>>::iterator ins_pos;
         list<PairedHit> chain;
         chain.push_back(move(open_hit));
         bool status;
         tie(ins_pos, status) = _open_mates.insert(make_pair(hit_id, chain));
         assert(status);
       }
       else{
         reset(orig_left, orig_right, orig_ref_id);
         LOG_ERR("POSSIBLE wrong alignment ", hit->ref_id(),":",hit->left()," with no gap between paired hits");
         return false;
       }
     } else{

       //if(iter_open->second.ref_id() != hit_ref_id)
         //return false;
       assert(!iter_open->second.empty());

       for(auto it = iter_open->second.begin(); it != iter_open->second.end(); ++it){
         bool strand_agree = it->strand() == hit_strand ||
              hit_strand == Strand_t::StrandUnknown ||
              it->strand() == Strand_t ::StrandUnknown;

         uint expected_pos = 0;
         if(it->_right_read){
            expected_pos = it->_right_read->partner_pos();
         }
         else{
            expected_pos = it->_left_read->partner_pos();
//#ifdef DEBUG
//     if(hit->left() == 46125){ //45791){
//       cout<<"it_strand"<<it->strand()<<endl;
//       cout<<"hit strand"<<hit_strand<<endl;
//
//     }
//#endif
         }
         if(it->left_pos() == hit_partner_pos &&
              it->ref_id() == hit_ref_id &&
              strand_agree &&
              expected_pos == hit->left())

         {
            if(it->_left_read == nullptr && it->_right_read){
              it->set_left_read(move(hit));
            }
            else if(it->_right_read == nullptr && it->_left_read){
             it->set_right_read(move(hit));
            }
            else{
              assert(false);
            }
            addHit(*it);
            iter_open->second.erase(it);
            if(iter_open->second.empty()){
              _open_mates.erase(iter_open);
            }
            return true;
         }

       }
       if(hit->partner_pos() > hit->right()){
         PairedHit open_hit(hit, nullptr);
         iter_open->second.push_back(open_hit);
       }
       else if(hit->partner_pos() < hit->left()){
         PairedHit open_hit(nullptr, hit);
         iter_open->second.push_back(open_hit);
       }
       else{
         reset(orig_left, orig_right, orig_ref_id);
         LOG_ERR("POSSIBLE wrong alignment ", hit->ref_id(),":",hit->left()," with no gap between paired hits");
         return false;
       }

     }
   }
   return true;
}

void HitCluster::clearOpenMates()
{
   _open_mates.clear();
}

int HitCluster::collapseHits()
{
   if(!_uniq_hits.empty()){
   //   cout<<this->left()<<"-"<<this->right()<<endl;
   //  cout<<_uniq_hits.size()<<endl;
     assert(false);
   }

   for (const auto& mates : _open_mates) {
      for (auto const& i : mates.second) _hits.push_back(i);
   }
   if(_hits.empty()) return 0;

   sort(_hits.begin(), _hits.end());

   for(size_t i = 0; i < _hits.size(); ++i){
     if(_uniq_hits.empty()){
       _uniq_hits.push_back(_hits[i]);
       _uniq_hits.back().add_2_collapse_mass(_hits[i].weighted_mass());
     }
     else{
       if(_uniq_hits.back() != _hits[i]){
         _uniq_hits.push_back(_hits[i]);
         _uniq_hits.back().add_2_collapse_mass(_hits[i].weighted_mass());
       }
       else{
         _uniq_hits.back().add_2_collapse_mass(_hits[i].weighted_mass());
       }
     }
   }

   return _uniq_hits.size();
}

void HitCluster::reweight_read(bool weight_bias)
{
   // this function is a place holder and does not
   // weight the reads
   assert(weight_bias == false);
   for(auto & hit:_hits){
     hit.init_raw_mass(); // set hit mass
     _weighted_mass += hit.raw_mass(); // set cluster mass
   }
   return;
}

//void HitCluster::reweight_read(const unordered_map<std::string, double>& kmer_bias, int num_kmers){
//
//   for(auto & hit:_hits){
//     hit.set_kmers(num_kmers);
//     hit.init_raw_mass();
//   }
//   for(auto & hit:_hits){
//     double weight = 0.0;
//     int num = 0;
//     for(auto const & kmer : hit._left_kmers){
//       auto got = kmer_bias.find(kmer.toString());
//       if(got != kmer_bias.end()){
//         weight += got->second;
//         num ++ ;
//       }
//     }
//
//     for(auto const & kmer : hit._right_kmers){
//       auto got = kmer_bias.find(kmer.toString());
//       if(got != kmer_bias.end()){
//         weight += got->second;
//         num ++ ;
//       }
//     }
//
//     double new_mass = hit.weighted_mass() * weight / num;
//     hit.weighted_mass(new_mass);
//     cout<<"raw mass "<<hit.raw_mass()<<endl;
//     cout<<"weighted mass "<<hit.weighted_mass()<<" and weight "<<weight<<endl;
//     if(weight == 0.0) exit(0);
//   }
//
//   for(auto const & hit:_hits){
//     _weighted_mass += hit.weighted_mass(); // set cluster mass
//   }
//
//}

double HitCluster::weighted_mass() const
{
   //cout<<_weighted_mass<<endl;
   return _weighted_mass;
}

void HitCluster::addWeightedMass(double m){
   _weighted_mass += m;
}

//void HitCluster::setBoundaries(){
///*
// * Call this after collapseHits
// */
//   if(enforce_ref_models && hasRefmRNAs()){
//     _leftmost = INT_MAX;
//     _rightmost = 0;
//     for(auto r : _ref_mRNAs){
//       _leftmost = min(_leftmost, r->left());
//       _rightmost = max(_rightmost, r->right());
//     }
//   }
//}


bool HitCluster::overlaps( const HitCluster& rhs) const{
   if(ref_id() != rhs.ref_id()) return false;
   return left() < rhs.left() ? right() >= rhs.left() : left() <= rhs.right();
}

//void HitCluster::guessStrand(){
//   int plus_intron_size =0;
//   int minus_intron_size =0;
//
//    for(auto r = _uniq_hits.cbegin(); r< _uniq_hits.cend(); ++r){
//      if(r->contains_splice()){
//         if(r->strand() == Strand_t::StrandPlus)
//           ++plus_intron_size;
//         else
//           ++minus_intron_size;
//      }
//    }
//    if(plus_intron_size == 0 && minus_intron_size >0){
//      _strand = Strand_t::StrandMinus;
//    }
//    else if(plus_intron_size > 0 && minus_intron_size ==0){
//      _strand = Strand_t::StrandPlus;
//    }
//    else if(plus_intron_size == 0 && minus_intron_size == 0){
//      _strand = Strand_t::StrandUnknown;
//    }
//    else{
//      if(plus_intron_size / minus_intron_size > _kMinFold4BothStrand)
//         _strand = Strand_t::StrandPlus;
//      if(minus_intron_size / plus_intron_size > _kMinFold4BothStrand)
//         _strand = Strand_t::StrandMinus;
//      _strand = Strand_t::StrandBoth;
//    }
//}

bool HitCluster::see_both_strands(){
//   bool see_plus = false;
//   bool see_minus = false;
//   for(auto &hit: _hits){
//     if(hit.contains_splice()){
//       if(hit.strand() == Strand_t::StrandPlus)
//         see_plus = true;
//       else
//         see_minus = true;
//     }
//   }
//   return see_plus && see_minus;
//   if(left() == 9092004){
//     cout<<"plus read: "<<_plus_strand_num_hits<<endl;
//     cout<<"minus read: "<<_minus_strand_num_hits<<endl;
//     exit(0);
//   }

   if(_plus_strand_num_hits >= 5 && _minus_strand_num_hits >= 5)
     return true;
   if(_plus_strand_num_hits == 0) return false;
   if(_minus_strand_num_hits == 0) return false;
   if(_minus_strand_num_hits < _plus_strand_num_hits &&
     _minus_strand_num_hits/(double) _plus_strand_num_hits > kMinIsoformFrac)
     return true;
   if(_plus_strand_num_hits < _minus_strand_num_hits &&
     _plus_strand_num_hits/(double) _minus_strand_num_hits > kMinIsoformFrac)
     return true;
   return false;
}

//int Sample::max_inner_dist() const
//{
//   if(_is_inspecting)
//     return kMaxInnerDist;
//   else{
//     if(_insert_size_dist->empty())
//       return kMaxInnerDist;
//     else
//       return _insert_size_dist->_end_offset -_hit_factory->_reads_table._read_len_abs*2;
//   }
//}
bool Sample::load_chrom_fasta(RefID seq_id)
{
   /*
   * Loading reference FASTA sequence without gene model
   */
   if (_fasta_getter == nullptr) throw std::runtime_error("ref fasta file not provided. Did \
                                                          you forget -f option?");
   string seq_name = _hit_factory->_ref_table.ref_real_name(seq_id);
   _fasta_interface->load2FaSeqGetter(*_fasta_getter, seq_name);
   return _fasta_getter->loadSeq();
}


string Sample::get_iso_seq(const shared_ptr<FaSeqGetter> &fa_getter, const Contig iso) const
{
   string iso_seq;
   for(size_t i = 0; i<iso._genomic_feats.size(); ++i){
     if(iso._genomic_feats[i]._match_op._code == S_MATCH){
       uint s = iso._genomic_feats[i].left();
       int l = iso._genomic_feats[i].len();
       iso_seq += fa_getter->fetchSeq(s, l);
     }
   }
   return iso_seq;
}



bool Sample::loadRefmRNAs(vector<unique_ptr<GffTree>> &gseqs, RefSeqTable &rt)
{
   cout<<"Has loaded transcripts from "<<gseqs.size()<<" Chromosomes/Scaffolds"<<endl;
   /*
   * Parse transcripts in GffTree to a vector of Contig objects.
   */

   //sort gseqs accroding to the observation order in ref_table
   // or if ref_table is empty, initialize it according to gseqs
   //ref_id in ref_table start from 0.
   if(rt.size() == 0){
     for(uint i=0; i<gseqs.size(); ++i){
       rt.get_id(gseqs[i]->_g_seq_name);
     }
   } else{
     for(uint i = 0; i<gseqs.size(); ++i){
       int idx = gseqs[i]->get_gseq_id();
       int ref_table_id = rt.get_id(gseqs[i]->_g_seq_name);
       if( idx != ref_table_id ){
         LOG_WARN("Warning: Sam file and Gff file are not sorted in the same order!");
         gseqs[i]->set_gseq_id(ref_table_id);
       }
     }
     sort(gseqs.begin(),gseqs.end(),
         [](const unique_ptr<GffTree> &lhs, const unique_ptr<GffTree> &rhs){
            return lhs->get_gseq_id() < rhs->get_gseq_id();
     });
   }
//   FaInterface fa_api;
//   FaSeqGetter *fsg = NULL;

// load sequencing if required
//   if(seqFile != NULL){
//     fa_api.initiate(seqFile);
//   }
   for(uint i = 0; i<gseqs.size(); ++i){// for loop for each chromosome
     GffTree * gseq = &(*gseqs[i]);
     int f = 0;
     int r = 0;
     int u = 0;
     RefID ref_id = rt.get_id(gseqs[i]->_g_seq_name);
     GffmRNA *mrna = NULL;
     vector<Contig> ref_mrna_for_chr;
//     if(fa_api.hasLoad()){
//       delete fsg;
//       fsg = NULL;
//       fsg = new FaSeqGetter();
//       fa_api.load2FaSeqGetter(*fsg,gseqs[i]->_g_seq_name);
//       if(fsg == NULL){
//         cerr<<"Reference sequence "<<gseqs[i]->_g_seq_name<<" can not be load!"<<endl;
//         cerr<<"Please check if the names of sequences in fasta file match the names in BAM/SAM file"<<endl;
//       }
//     }
     int f_total = gseqs[i]->_forward_rnas.size();
     int r_total = gseqs[i]->_reverse_rnas.size();
     int u_total = gseqs[i]->_unstranded_rnas.size();
     while(!(f==f_total && r==r_total && u == u_total)){
       Strand_t strand;
       if(f <f_total){
         mrna = &(*gseqs[i]->_forward_rnas[f++]);
         strand = Strand_t::StrandPlus;
       } else if(r < r_total){
         mrna = &(*gseqs[i]->_reverse_rnas[r++]);
         strand = Strand_t::StrandMinus;
       } else{
         mrna = &(*gseqs[i]->_unstranded_rnas[u++]);
         strand = Strand_t::StrandUnknown;
       }
       if(mrna->_exons.size() == 0){
         continue;
       }
       vector<GenomicFeature> feats;
       for(uint e = 0; e < mrna->_exons.size(); ++e){
         GffExon& ex = *(mrna->_exons[e]);
         feats.push_back(GenomicFeature(Match_t::S_MATCH, ex._iv.left(), ex._iv.right()-ex._iv.left()+1));
         if( e + 1 < mrna->_exons.size()){
            GffExon& next_ex  = *(mrna->_exons[e+1]);
            feats.push_back(GenomicFeature(Match_t::S_INTRON, ex._iv.right()+1, next_ex._iv.left()-1-ex._iv.right() ));
         }
       }
       Contig ref_contig(ref_id, 0, strand,1.0, feats, true);
       ref_contig.annotated_trans_id(mrna->_transcript_id);
       ref_contig.parent_id() = mrna->getParentGene()->_gene_id;
       ref_contig.mass(1.0);
       //cout<<"ref contig left pos "<<ref_contig.left()<<endl;
       ref_mrna_for_chr.push_back(ref_contig);
     }// end while loop
     sort(ref_mrna_for_chr.begin(), ref_mrna_for_chr.end());
     _ref_mRNAs.insert(_ref_mRNAs.end(), ref_mrna_for_chr.begin(), ref_mrna_for_chr.end());

     ref_mrna_for_chr.clear();
   }//end for loop

   //delete fsg;
   //fsg = NULL;
   return true;
}

double Sample::next_valid_alignment(ReadHit& readin){
   const char* hit_buf=NULL;
   size_t hit_buf_size = 0;
   double raw_mass = 0.0;
   while(true){
     if(!_hit_factory->nextRecord(hit_buf, hit_buf_size)) {
       break;
     }

     if(!_hit_factory->getHitFromBuf(hit_buf, readin)) continue;
     if(readin.ref_id() == -1) continue; // unmapped read
     raw_mass += readin.mass(); // suck in read mass for future if mask_gtf is used.

     if(_prev_hit_ref_id != -1){
       if(_prev_hit_ref_id > readin.ref_id() ||
            ( _prev_hit_ref_id == readin.ref_id() && _prev_hit_pos > readin.left())
         )
       {
         const string cur_chr_name = _hit_factory->_ref_table.ref_real_name(readin.ref_id());
         const string last_chr_name = _hit_factory->_ref_table.ref_real_name(_prev_hit_ref_id);
         if(_is_inspecting){
            LOG_ERR("BAM file not sort correctly!");
            LOG_ERR("The current position is: ", cur_chr_name, ":", readin.left());
            LOG_ERR("and previous position is: ", last_chr_name, ":", _prev_hit_pos);
         }
       }
     }

     _prev_hit_ref_id = readin.ref_id();
     _prev_hit_pos = readin.left();
     break;
   }
   _hit_factory->_reads_table._read_len_abs = max(_hit_factory->_reads_table._read_len_abs, readin.read_len());
   return raw_mass;
}

double Sample::rewindHit(const ReadHit& rh)
{
   double mass = rh.mass();
   _hit_factory->undo_hit();
   return mass;
}

int Sample::addRef2Cluster(HitCluster &cluster_out){
   if(_refmRNA_offset >=  _ref_mRNAs.size()) {
     _has_load_all_refs = true;
     return 0;
   }

   // add first rna
   //cout<<"offset: "<<_refmRNA_offset<<endl;
   //cout<<"_ref_mRNAs size: "<<_ref_mRNAs[0].parent_id()<<endl;

   cluster_out.gene_id() = _ref_mRNAs[_refmRNA_offset].parent_id();
   cluster_out.addRefContig(&_ref_mRNAs[_refmRNA_offset++]);
   if(_refmRNA_offset >= _ref_mRNAs.size()){
     _has_load_all_refs = true;
     return 1;
   }

   // add the rest if overlapped with first
   size_t i = 0;
   while(i < cluster_out._ref_mRNAs.size()){
     Contig* ref = cluster_out._ref_mRNAs[i];
     if (Contig::overlaps_directional(*ref, _ref_mRNAs[_refmRNA_offset])) {
       cluster_out.addRefContig(&_ref_mRNAs[_refmRNA_offset++]);
       if(_refmRNA_offset >= _ref_mRNAs.size()){
         _has_load_all_refs = true;
         return cluster_out._ref_mRNAs.size();
       }
       i=0;
     }
     else{
       ++i;
     }
   }
   return cluster_out._ref_mRNAs.size();
}

void Sample::rewindReference(HitCluster &clusterOut , int num_regress)
{
   clusterOut.left(UINT_MAX);
   clusterOut.right(0);
   clusterOut.ref_id(-1);
   clusterOut._ref_mRNAs.clear();
   _refmRNA_offset -= num_regress;
   assert(_refmRNA_offset >= 0);
}

void Sample::reset_refmRNAs()
{
   _refmRNA_offset = 0;
   _has_load_all_refs = false;
   if (!no_assembly) {
     _ref_mRNAs.clear();
     move(_assembly.begin(), _assembly.end(), back_inserter(_ref_mRNAs));
     _assembly.clear();
     sort(_ref_mRNAs.begin(), _ref_mRNAs.end());
   }
}

int Sample::nextCluster_denovo(HitCluster &clusterOut,
                              uint next_ref_start_pos,
                              RefID next_ref_start_ref
                              )
{
   if(!_hit_factory->recordsRemain()) return -1;
   while(true){
     ReadHitPtr new_hit(new ReadHit());
     double mass = next_valid_alignment(*new_hit);
     //cout<<new_hit->left()<<endl;

     if(!_hit_factory->recordsRemain()){
       return clusterOut.size();
     }


     if(new_hit->ref_id() > next_ref_start_ref ||
      (new_hit->ref_id() == next_ref_start_ref && new_hit->right() >= next_ref_start_pos)){
       rewindHit(*new_hit);
       return clusterOut.size();
     }

     if(clusterOut.ref_id() == -1){ // add first hit

       clusterOut.addOpenHit(new_hit, true, true);
       clusterOut.addRawMass(mass);
     } else { //add the rest
       if(hit_lt_cluster(*new_hit, clusterOut, kMaxOlapDist)){
         // should never reach here
         LOG_ERR("In alignments.cpp: It appears that SAM/BAM not sorted!");
       }
       if(hit_gt_cluster(*new_hit, clusterOut, kMaxOlapDist)){
         // read has gone to far.
         rewindHit(*new_hit);
         break;
       }
       clusterOut.addOpenHit(new_hit, true, true);
       clusterOut.addRawMass(mass);
     }
   }
   return clusterOut.size();
}

int Sample::nextClusterRefDemand(HitCluster &clusterOut){
   int wiggle_room = utilize_ref_models ? 0 : kMaxOlapDist;

   if (!hasLoadRefmRNAs()) {
     std::cerr<<"if you use --no-assembly option, you must provide gff file through -g option!"<<std::endl;
     assert(false);
   }
   if (!_hit_factory->recordsRemain()) return -1;
   int num_added_refmRNA = addRef2Cluster(clusterOut);

   if (num_added_refmRNA == 0) {
     return -1;
   }

   //cout<<"cluster: "<<clusterOut.left()<<"-"<<clusterOut.right()<<endl;
   while (true) {
     if(!_hit_factory->recordsRemain()){
       break;
     }
     ReadHitPtr new_hit(new ReadHit());
     double mass = next_valid_alignment(*new_hit);
     if (hit_lt_cluster(*new_hit, clusterOut, wiggle_room)) {  //hit hasn't read this region

     } else if (hit_gt_cluster(*new_hit, clusterOut, wiggle_room)) {
       rewindHit(*new_hit);
       break;
     } else {
       //cout<<"new_hit :"<<new_hit->left()<<"-"<<new_hit->right()<<endl;
       clusterOut.addOpenHit(new_hit, false, false);
       clusterOut.addRawMass(mass);
     }
   }  //end while loop
   return clusterOut.size();
}

void Sample::preProcess(FILE *log) {
   _is_inspecting = true;
   curr_thread_num = 0;
   const RefSeqTable & ref_t = _hit_factory->_ref_table;

   _num_cluster = 0;
   while(true){
     shared_ptr<HitCluster> cluster (new HitCluster());
     if(-1 == nextClusterRefDemand(*cluster)){
       break;
     }
     if(cluster->ref_id() == -1) continue;
     cluster->_id = ++_num_cluster;
#if ENABLE_THREADS
     if(use_threads){
       while(true){
         if(curr_thread_num < num_threads){
            break;
         }
         this_thread::sleep_for(chrono::milliseconds(3));
       }
       ++curr_thread_num;
       thread worker ([=] {
         finalizeCluster(cluster, true);
         fragLenDist(ref_t, cluster->ref_mRNAs(), cluster, log);
         --curr_thread_num;
       });
       worker.detach();
     }else {
       finalizeCluster(cluster, true);
       fragLenDist(ref_t, cluster->ref_mRNAs(), cluster, log);
     }
#endif
   } //end while(true)

#if ENABLE_THREADS
   if(use_threads){
     while(true){
       if(curr_thread_num==0){
         break;
       }
       this_thread::sleep_for(chrono::milliseconds(5));
     }
   }
#endif
}

int Sample::nextCluster_refGuide(HitCluster &clusterOut)
{
   //bool skip_read = false;
   if(!_hit_factory->recordsRemain()) return -1;
   int wiggle_room = utilize_ref_models ? 0 : kMaxOlapDist;

   if(!hasLoadRefmRNAs()){
     return nextCluster_denovo(clusterOut);
   }
   else{
     int num_added_refmRNA = addRef2Cluster(clusterOut);
     //if all ref mRNAs have been loaded but alignments haven't
     if( num_added_refmRNA == 0){
       return nextCluster_denovo(clusterOut);
     }
     //else add as many as alignment possible
     while(true){
       ReadHitPtr new_hit(new ReadHit());
       double mass = next_valid_alignment(*new_hit);

       if(hit_lt_cluster(*new_hit, clusterOut, wiggle_room)){ // hit hasn't reach reference region
         rewindHit(*new_hit);
         if(_has_load_all_refs){
            rewindReference(clusterOut, num_added_refmRNA);
            return nextCluster_denovo(clusterOut);
         } else{
#ifdef DEBUG
            assert(_ref_mRNAs[_refmRNA_offset].featSize() > 0);
#endif
            uint next_ref_start_pos = _ref_mRNAs[_refmRNA_offset].left();
            uint next_ref_start_ref = _ref_mRNAs[_refmRNA_offset].ref_id();
   //#ifdef DEBUG
   //         cout<<"next:"<<next_ref_start_pos<<"\t hit chr:"<<new_hit->ref_id()<<"\thit left "<<new_hit->left()<<"\t previous:\t"<<clusterOut.left()<<endl;
   //#endif
            rewindReference(clusterOut, num_added_refmRNA);
            return nextCluster_denovo(clusterOut, next_ref_start_pos, next_ref_start_ref);
         }
       }

       if(hit_gt_cluster(*new_hit, clusterOut, wiggle_room)){ // read has gone too far.
         rewindHit(*new_hit);
         break;
       }

       clusterOut.addOpenHit(new_hit, false, false);
       clusterOut.addRawMass(mass);
       if (!_hit_factory->recordsRemain()) return clusterOut.size();
     } // end while loop
   } // end loadRefmRNAs
   return clusterOut.size();
}

void Sample::reAssignClusters(HitCluster & last, HitCluster &cur){
   /*
   * reassign paired hits;
   */
   sort(last._hits.begin(), last._hits.end());
   vector<PairedHit> imcompatible_hits;
   vector<bool> is_imcomp(last._hits.size(), false);
   bool first_incompatible = false;
   for(uint i=0; i != last._hits.size();++i){
     PairedHit* hit_it = &last._hits[i];
     if(hit_it->contains_splice()){
       if(hit_it->strand() != last.guessStrand()){
         first_incompatible = true;
         imcompatible_hits.push_back(*hit_it);
         is_imcomp[i] = true;
       }
       else{
         ++hit_it;
         first_incompatible = false;
       }
     }
     else{
       if(first_incompatible){
         imcompatible_hits.push_back(*hit_it);
         is_imcomp[i] = true;
       }
       else{
         ++hit_it;
       }
     }
   }
   for(auto it = imcompatible_hits.cbegin(); it != imcompatible_hits.cend(); ++it){
     if(hit_complete_within_cluster(*it, cur, 0))
       cur.addHit(*it);
   }

   vector<PairedHit> reduced_hits;
   for(uint i=0; i != last._hits.size(); ++i){
     //IF DEBUG
     //cout<<last.left()<<":"<<last._hits.back().left_pos()<<endl;
     //ENDIF
     if(is_imcomp[i] == false) reduced_hits.push_back(last._hits[i]);
   }
   last._hits = reduced_hits;

   /*
   * reassign open hits;
   */
   for(auto it = last._open_mates.begin(); it != last._open_mates.end(); ++it){
     for(auto hit = it->second.begin(); hit != it->second.end(); ++hit){
       if(hit->_left_read){
         cur.addOpenHit(hit->_left_read,false,false);
       }
       else{
         assert(hit->_right_read);
         cur.addOpenHit(hit->_right_read,false,false);
       }
     }

   }
}


void Sample::finalizeCluster(shared_ptr<HitCluster> cluster, bool clear_open_mates ){
   if (cluster->size() == 0) {
      return;
   }
   if(clear_open_mates){
     cluster->clearOpenMates();
   }
   cluster->reweight_read(false);
   cluster->collapseHits();
   //cluster->setBoundaries(); // set boundaries if reference exist.
}

void Sample::fragLenDist(const RefSeqTable &ref_t,
               const std::vector<Contig> &transcripts,
               const shared_ptr<HitCluster> cluster,
               FILE *plogfile) {

   if (transcripts.empty()) return;
   _total_mapped_reads += (int) cluster->weighted_mass();

   vector<Contig> hits;
   for (auto r = cluster->_uniq_hits.cbegin(); r != cluster->_uniq_hits.cend(); ++r) {
     Contig hit(*r);
     hits.push_back(hit);
   }

   if (transcripts.size() == 1) {
     for (auto const &assembled_transcript: transcripts) {
       for (size_t h = 0; h < hits.size(); ++h) {
         if (hits[h].is_single_read()) continue;
         if (!Contig::is_compatible(hits[h], assembled_transcript)) continue;
         double frag_len = Contig::exonic_overlaps_len(assembled_transcript,
                                    hits[h].left(),
                                    hits[h].right());

#if ENABLE_THREADS
         if (use_threads) {
            thread_pool_lock.lock();
            _hit_factory->_reads_table._frag_dist.push_back(frag_len);
            thread_pool_lock.unlock();
         } else {
            _hit_factory->_reads_table._frag_dist.push_back(frag_len);
         }
#else
         _hit_factory->_reads_table._frag_dist.push_back(frag_len);
#endif
       }// end for
     }// end for
   } // end if (assembled_transcripts.size()==1)

#if ENABLE_THREADS
   if (use_threads) {
     out_file_lock.lock();
   }
#endif
   fprintf(plogfile, "Finish inspecting locus: %s:%d-%d\n", ref_t.ref_real_name(cluster->ref_id()).c_str(),
       cluster->left(), cluster->right());
   fprintf(plogfile, "Found %d of ref mRNAs from the reference gtf file.\n", cluster->_ref_mRNAs.size());
   fprintf(plogfile, "Number of total unique hits: %d\n\n", cluster->_uniq_hits.size());

#if ENABLE_THREADS
   if (use_threads) {
     out_file_lock.unlock();
   }
#endif
}

vector<Contig> Sample::assembleCluster(const RefSeqTable &ref_t, shared_ptr<HitCluster> cluster, FILE *plogfile) {
   /*
    * Isoform id is integer from 1. The index
    * for each isoform in vector<Isoform> is id+1.
    * */


   //local variables
   vector<Contig> assembled_transcripts;
   vector<float> exon_doc;
   IntronMap intron_counter;
   vector<GenomicFeature> exons;
   vector<Contig> hits;

   if (cluster->size() == 0) {
//#if ENABLE_THREADS
//      if (use_threads) decr_pool_count();
//#endif
      return assembled_transcripts;
   }

   // if single-end library not need to calculate fragment length distribution
   if (SINGLE_END_EXP && _is_inspecting) {
//#if ENABLE_THREADS
//      if (use_threads)
//         decr_pool_count();
//#endif
      return assembled_transcripts;
   }

   //enough read for calculating empirical distribution of reads.
   if (_is_inspecting &&
      _hit_factory->_reads_table._frag_dist.size() > kMaxReadNum4FD) {
//#if ENABLE_THREADS
//      if (use_threads)
//         decr_pool_count();
//#endif
      return assembled_transcripts;
   }


   if (cluster->len() > kMaxGeneLength) {
//#if ENABLE_THREADS
//      if (use_threads)
//         decr_pool_count();
//#endif
      return assembled_transcripts;
   }

   cluster->_strand = cluster->guessStrand();


   for (auto r = cluster->_uniq_hits.cbegin(); r != cluster->_uniq_hits.cend(); ++r) {
      Contig hit(*r);
      hits.push_back(hit);
   }

   if (cluster->hasRefmRNAs() && utilize_ref_models) {
      for (auto i: cluster->_ref_mRNAs)
         hits.push_back(*i);
      hits.back()._is_ref = true;
      //  avg_dep = compute_doc(cluster->left(), cluster->right(), hits, exon_doc, intron_counter, kMaxSmallAnchor);
   }
   sort(hits.begin(), hits.end());
   size_t s = cluster->right() - cluster->left() + 1;
   exon_doc.resize(s, 0);
   double avg_dep = 0.0;
   avg_dep = compute_doc(cluster->left(), cluster->right(), hits, exon_doc, intron_counter, kMaxSmallAnchor);
   //cout<<avg_dep<<endl;

   if (avg_dep < kMinDepth4Locus) {
//#if ENABLE_THREADS
//      if (use_threads)
//         decr_pool_count();
//#endif
      return assembled_transcripts;
   }

   filter_intron(cluster->left(), exon_doc, intron_counter);

   //local variables
   FlowNetwork flow_network;
   Graph::NodeMap<const GenomicFeature *> node_map(flow_network._g);
   Graph::ArcMap<int> cost_map(flow_network._g);
   Graph::ArcMap<int> min_flow_map(flow_network._g);
   vector<vector<Graph::Arc>> path_cstrs;
   vector<vector<GenomicFeature>> assembled_feats;
   vector<vector<size_t>> constraints;

//#ifdef DEBUG
//     cout<<"cluster starts at: "<<cluster->left()<<endl;
//     cout<<"exon coverage:"<<endl;
//     for(auto i : exon_doc)
//       cout<<i;
//     cout<<endl;
//#endif

   flow_network.splicingGraph(cluster->ref_id(), cluster->left(), exon_doc, intron_counter, exons);

//#ifdef DEBUG
//              cout<<"frag_len: "<<frag_len<<endl;
//              cout<<"frag location>> "<<hits[h].ref_id()<<":"<<hits[h].left()<<"-"<<hits[h].right()<<endl;
//#endif

   constraints = flow_network.findConstraints(exons, hits);

   bool stat = flow_network.createNetwork(hits, exons, intron_counter,
                                 constraints, node_map, cost_map, min_flow_map, path_cstrs);
   if (!stat) {
//#if ENABLE_THREADS
//      if (use_threads)
//         decr_pool_count();
//#endif
      return assembled_transcripts;
   }

   bool stat2 = flow_network.solveNetwork(node_map, exons, path_cstrs, cost_map, min_flow_map, assembled_feats);
   if (!stat2) {
//#if ENABLE_THREADS
//      if (use_threads)
//         decr_pool_count();
//#endif
      return assembled_transcripts;
   }

   assembled_transcripts = assemble_2_contigs(assembled_feats, cluster->ref_id(), cluster->strand());

   this->fragLenDist(ref_t, assembled_transcripts, cluster, plogfile);

   /* Merge transfrags if they
    *  do not overlaps.
    *  Begin*/
   if (kCombineShrotTransfrag && assembled_transcripts.size() > 1) {
      vector<pair<uint, uint>> to_merge;
      vector<Contig> merged_transcripts;
      for (uint i = 0; i < assembled_transcripts.size() - 1; ++i) {
         for (uint j = i + 1; j < assembled_transcripts.size(); ++j) {
            if (!Contig::overlaps_directional(assembled_transcripts[i], assembled_transcripts[j])) {
               to_merge.emplace_back(i, j);
            }
         }
      }
      if (!to_merge.empty()) {
         vector<bool> flags(assembled_transcripts.size(), false);
         for (auto iso_pair: to_merge) {
            uint i = iso_pair.first;
            uint j = iso_pair.second;
            flags[i] = true;
            flags[j] = true;
            Contig new_contig;
            if (assembled_transcripts[i] < assembled_transcripts[j]) {
               new_contig = assembled_transcripts[i];
               GenomicFeature intron(Match_t::S_INTRON, new_contig.right() + 1, \
                assembled_transcripts[j].left() - new_contig.right() - 1);
               new_contig._genomic_feats.push_back(intron);
               new_contig._genomic_feats.insert(new_contig._genomic_feats.end(), \
                                    assembled_transcripts[j]._genomic_feats.begin(), \
                                    assembled_transcripts[j]._genomic_feats.end());
            } else {
               new_contig = assembled_transcripts[j];
               GenomicFeature intron(Match_t::S_INTRON, new_contig.right() + 1, \
                assembled_transcripts[i].left() - new_contig.right() - 1);
               new_contig._genomic_feats.push_back(intron);
               new_contig._genomic_feats.insert(new_contig._genomic_feats.end(), \
                                    assembled_transcripts[i]._genomic_feats.begin(), \
                                    assembled_transcripts[i]._genomic_feats.end());
            }
            merged_transcripts.push_back(new_contig);
         } // end for

         for (uint i = 0; i < assembled_transcripts.size(); ++i) {
            if (flags[i]) continue;
            merged_transcripts.push_back(assembled_transcripts[i]);
         }
         assembled_transcripts = merged_transcripts;
      } // end if(!to_merge.empty())
   }  //end if(kCombineShrotTransfrag)
   /* Merge transfrags if they
    *  do not overlaps.
    *  End*/
   return assembled_transcripts;
}

void Sample::quantifyCluster(const RefSeqTable &ref_t, const shared_ptr<HitCluster> cluster,
                 const vector<Contig> &assembled_transcripts, FILE *pfile, FILE *plogfile) const {

   vector<Isoform> isoforms;
   map<int, int> iso_2_len_map;

   vector<GenomicFeature> exons; //= Contig::uniqueFeatsFromContigs(assembled_transcripts, Match_t::S_MATCH);
   vector<Contig> hits;

   for (auto r = cluster->_uniq_hits.cbegin(); r != cluster->_uniq_hits.cend(); ++r) {
     Contig hit(*r);
     hits.push_back(hit);
   }
   //cout<<"open hit "<<cluster->_open_mates.size()<<endl;

   assert(assembled_transcripts.size());
   // prepare exon segments
   for(const auto &t: assembled_transcripts) {
     for (const auto &f: t._genomic_feats) {
       if (f._match_op._code == Match_t::S_MATCH) exons.push_back(f);
     }
   }
   sort(exons.begin(), exons.end());
   auto last = unique(exons.begin(), exons.end());
   exons.erase(last, exons.end());
   IRanges<GenomicFeature, false> exons_iranges(exons);
   vector<GenomicFeature> reduced_exons = exons_iranges.disjoint();
   exons = reduced_exons;

   for(const auto &t: assembled_transcripts){
     Isoform iso(exons, t, t.parent_id(), t.annotated_trans_id(), cluster->_id);
     int idx = PushAndReturnIdx<Isoform>(iso, isoforms);
     iso_2_len_map[idx] = t.exonic_length();
   }
   LocusContext est(_insert_size_dist, _hit_factory->_reads_table._read_len_abs,
                    plogfile, hits, isoforms, exons, _fasta_getter);

   for (size_t i = 0; i < est.num_hit(); ++i) {
      for (size_t k = 0; k < est.num_isoform(); ++k) {
         for (size_t l = 0; l < est.num_exon_bin(); ++l) {
            cout<<est.GetXikl(i,k,l)<<endl;
         }
      }
   }

   exit(0);
   //est.set_empirical_bin_weight(iso_2_bins_map, iso_2_len_map, cluster->collapse_mass(), exon_bin_map);
   est.set_theory_bin_weight(iso_2_len_map, isoforms);
   //est.calculate_raw_iso_counts(iso_2_bins_map, exon_bin_map);
   bool success = est.estimate_abundances(this->total_mapped_reads(), \
                                iso_2_len_map, isoforms, BIAS_CORRECTION);
   if(success){
#if ENABLE_THREADS
     if(use_threads)
       out_file_lock.lock();
#endif
     /* Print locus coordinates*/
     cerr<<ref_t.ref_real_name(cluster->ref_id())<<"\t"<<cluster->left()<<"\t"<<cluster->right()<<" finishes abundances estimation"<<endl;
     for(auto & iso: isoforms){
       iso._contig.print2gtf(pfile, _hit_factory->_ref_table, iso._FPKM_s,
                    iso._TPM_s, iso._gene_str, iso._isoform_str);
     }
     fprintf(plogfile, "Finish abundances estimation at locus: %s:%d-%d\n", ref_t.ref_real_name(cluster->ref_id()).c_str(), cluster->left(), cluster->right());
   }
#if ENABLE_THREADS
   if(use_threads) {
     out_file_lock.unlock();
     decr_pool_count();
   }
#endif
}

void Sample::addAssembly(std::vector<Contig>& assembs, int cluster_id) {
   int tid=0;
   for (Contig& asmb: assembs) {
     ++tid;
     asmb.parent_id() = "gene."+to_string(cluster_id);
     asmb.annotated_trans_id("transcript." + to_string(cluster_id) + "." + to_string(tid));
   }

#if ENABLE_THREADS
   if (use_threads) thread_pool_lock.lock();
#endif

   if (!assembs.empty()) std::move(assembs.begin(), assembs.end(), std::back_inserter(_assembly));

#if ENABLE_THREADS
   if(use_threads) thread_pool_lock.unlock();
#endif
}

void Sample::inspectSample(FILE *plogfile)
/*
 *  First run-through to calculate fragment distribution FD
 * */
{
   _is_inspecting = true;
   const RefSeqTable & ref_t = _hit_factory->_ref_table;
   shared_ptr<HitCluster> last_cluster (new HitCluster());
   if( -1 == nextCluster_refGuide(*last_cluster) ) {
     return;
   }
   _num_cluster = 1;

   RefID current_ref_id = last_cluster->ref_id();
   if(BIAS_CORRECTION)
     load_chrom_fasta(current_ref_id);
   _current_chrom = ref_t.ref_real_name(last_cluster->ref_id());

   while(true){
     shared_ptr<HitCluster> cur_cluster (new HitCluster());
     if (!last_cluster->hasRefmRNAs() && last_cluster->see_both_strands()){
       cur_cluster->left(last_cluster->left());
       cur_cluster->right(last_cluster->right());
       cur_cluster->ref_id(last_cluster->ref_id());
       reAssignClusters(*last_cluster, *cur_cluster);
     } else{
       if (-1 == nextCluster_refGuide(*cur_cluster)){
         break;
       }
     }


     if (last_cluster->ref_id() == -1) {
       last_cluster = move(cur_cluster);
       continue;
     }

//Begin loading ref seqs
     if(current_ref_id != last_cluster->ref_id()){
       current_ref_id = last_cluster->ref_id();
       if(BIAS_CORRECTION){
#if ENABLE_THREADS
         if(use_threads){
            while(true){
              if(curr_thread_num==0){
                break;
              }
              this_thread::sleep_for(chrono::milliseconds(3));
            }
            load_chrom_fasta(current_ref_id);
         }
#endif
         load_chrom_fasta(current_ref_id);
       }
     }
     if(_current_chrom != ref_t.ref_real_name(last_cluster->ref_id())){
       _current_chrom = ref_t.ref_real_name(last_cluster->ref_id());
     }
//End loading ref seqs
     last_cluster->_id = ++_num_cluster;
#if ENABLE_THREADS
     if(use_threads){
       while(true){
         if(curr_thread_num < num_threads){
            break;
         }
         this_thread::sleep_for(chrono::milliseconds(3));
       }
       ++curr_thread_num;
       thread worker ([=] {
            finalizeCluster(last_cluster, true);
            vector<Contig> asmb = this-> assembleCluster(ref_t, last_cluster, plogfile);
            this->addAssembly(asmb, last_cluster->_id);
            --curr_thread_num;
            });
       //thread worker{&Sample::finalizeAndassembleCluster, this, ref_t, last_cluster, NULL, NULL};
       worker.detach();
     }else{
       finalizeCluster(last_cluster, true);
       vector<Contig> asmb = assembleCluster(ref_t, last_cluster, plogfile);
       this->addAssembly(asmb, last_cluster->_id);
     }
#else
   finalizeCluster(last_cluster, true);
   assembleCluster(ref_t, last_cluster, plogfile);
   this->addAssembly(asmb);
#endif
     // _total_mapped_reads += (int) last_cluster->raw_mass();
     //cout<<"weighted cluster mass: "<<last_cluster->_weighted_mass<<endl;
     fprintf(plogfile, "Inspect gene: %s:%d-%d\n", ref_t.ref_real_name(last_cluster->ref_id()).c_str(), last_cluster->left(), last_cluster->right());
     fprintf(plogfile, "Has inspected %d reads\n", (int)_total_mapped_reads);
     last_cluster = move(cur_cluster);
   }

//make sure all threads have finished
#if ENABLE_THREADS
   if(use_threads){
     while(true){
#ifdef DEBUG
       //cout<<curr_thread_num<<endl;
#endif
       if(curr_thread_num==0){
         break;
       }
       this_thread::sleep_for(chrono::milliseconds(3));
     }
   }
#endif
   last_cluster->_id = ++_num_cluster;
   finalizeCluster(last_cluster, true);
   vector<Contig> asmb = assembleCluster(ref_t, last_cluster, plogfile);
   this->addAssembly(asmb, last_cluster->_id);
   fprintf(plogfile, "Inspect gene: %s:%d-%d\n", ref_t.ref_real_name(last_cluster->ref_id()).c_str(), last_cluster->left(), last_cluster->right());
   fprintf(plogfile, "Has inspected %d reads\n", (int)_total_mapped_reads);
   //_total_mapped_reads += (int)last_cluster->raw_mass();
#if ENABLE_THREADS
   if(use_threads)
     curr_thread_num = 0;
#endif
   //_is_inspecting = false;
}

int Sample::total_mapped_reads() const
{
   return _total_mapped_reads;
}

void Sample::procSample(FILE *pfile, FILE *plogfile)
{
/*
 * The major function which calls nextCluster() and finalizes cluster and
 * assemble each cluster-> Right now only nextCluster_refGuide() is implemented.
 * if no reference mRNA than nextCluster_refGuide will call nextCluster_denovo()
 */
   //cout<<"reach in procSample"<<endl;
   _is_inspecting = false;
   _hit_factory->reset();
   reset_refmRNAs();
   const RefSeqTable & ref_t = _hit_factory->_ref_table;
   int current_ref_id = INT_MAX;

   while(true){
     //++_num_cluster;
     shared_ptr<HitCluster> cluster (new HitCluster());
     if(-1 == nextClusterRefDemand(*cluster)){
       break;
     }
     if(cluster->ref_id() == -1) continue;

     //cout<<"uniq "<<cluster->_uniq_hits.size()<<" cluster size "<<cluster->size()<<endl;
      //begin load fasta genome
     if(current_ref_id != cluster->ref_id()){
       current_ref_id = cluster->ref_id();
       if(BIAS_CORRECTION){
#if ENABLE_THREADS
         if(use_threads){
            while(true){
              if(curr_thread_num==0){
                break;
              }
              this_thread::sleep_for(chrono::milliseconds(3));
            }
            load_chrom_fasta(current_ref_id);
         }
#endif
         load_chrom_fasta(current_ref_id);
       }
     }
     //end load fasta genome

#if ENABLE_THREADS
     if(use_threads){
       while(true){
         if(curr_thread_num < num_threads){
            break;
         }
         this_thread::sleep_for(chrono::milliseconds(3));
       }
       ++curr_thread_num;
       thread worker ([=] {
            finalizeCluster(cluster, false);
            this->quantifyCluster(ref_t, cluster, cluster->ref_mRNAs(), pfile, plogfile);
       });
       worker.detach();
     }else {
       finalizeCluster(cluster, false);
       this->quantifyCluster(ref_t, cluster, cluster->ref_mRNAs(), pfile, plogfile);
     }
#else
     finalizeCluster(last_cluster, true);
     this->quantifyCluster(ref_t, last_cluster, transcripts, pfile, plogfile);
#endif
   } //end while(true)

#if ENABLE_THREADS
   if(use_threads){
     while(true){
       if(curr_thread_num==0){
         break;
       }
       this_thread::sleep_for(chrono::milliseconds(5));
     }
   }
#endif

}


double Sample::compute_doc(const uint left, const uint right, const vector<Contig> & hits,
     vector<float> &exon_doc,  IntronMap &intron_counter, uint smallOverHang)
{

   assert(right > left);
   for(size_t i = 0; i<hits.size(); ++i){
     const vector<GenomicFeature> & g_feats = hits[i]._genomic_feats;
     for(size_t j = 0; j<g_feats.size(); ++j){
       const GenomicFeature & gf = g_feats[j];
       if( gf._match_op._code == Match_t::S_MATCH){
         size_t l  = gf.left() < left ? left: gf.left();
         size_t r = gf.right() > right ? right: gf.right();
         for(size_t p = l; p < r+1; ++p){
            exon_doc[p-left] += hits[i].mass();
         }
       }
       else if( gf._match_op._code == Match_t::S_INTRON){
         if(gf.left() < left || gf.right() > right)
            continue;
         IntronTable cur_intron(gf.left(), gf.right());
         pair<uint,uint> coords(cur_intron.left, cur_intron.right);
         if(intron_counter.empty()){
            cur_intron.total_junc_reads += hits[i].mass();
            if(g_feats[j-1]._match_op._len < smallOverHang ||
                g_feats[j+1]._match_op._len < smallOverHang){
              cur_intron.small_span_read += hits[i].mass();
            }
            intron_counter.emplace(coords, cur_intron);
            continue;
         }
         auto it = intron_counter.find(coords);
         if( it != intron_counter.end()){
            it->second.total_junc_reads += hits[i].mass();
            if(g_feats[j-1]._match_op._len < smallOverHang ||
                g_feats[j+1]._match_op._len < smallOverHang){
              it->second.small_span_read += hits[i].mass();
            }
         }
         else{
            cur_intron.total_junc_reads += hits[i].mass();
            if(g_feats[j-1]._match_op._len < smallOverHang ||
                g_feats[j+1]._match_op._len < smallOverHang){
              cur_intron.small_span_read += hits[i].mass();
            }
            intron_counter.emplace(coords, cur_intron);
         }
       }
     }
   }
   int num_nt = 0;
   for(uint i=0; i != exon_doc.size(); ++i){
     if(exon_doc[i] > 0) ++num_nt;
   }
   if(num_nt == 0) return 0.0;
   double total_depth = accumulate(exon_doc.begin(), exon_doc.end(), 0.0);
   return total_depth/num_nt;
}

void Sample::filter_intron(uint cluster_left,
     vector<float> &exon_doc, IntronMap& intron_counter)
{
   // bad_intron_pos for indexes for intron to be drop off in intron_counter
   vector<float> intron_doc(exon_doc.size(),0.0);

//Filtering one: by overlapping with better intron
   vector<pair<uint, uint>> bad_intron_pos;
   for(auto i = intron_counter.cbegin(); i != intron_counter.cend(); ++i){
     for(auto j = next(intron_counter.cbegin()); j != intron_counter.cend(); ++j){
       if(IntronTable::overlap(i->second, j->second)){
          float depth_i = i->second.total_junc_reads;
          float depth_j = j->second.total_junc_reads;
         if( depth_i / depth_j < kMinIsoformFrac){
            bad_intron_pos.push_back(i->first);
            if(_is_inspecting){
              LOG("Filtering overlapping intron by depth: ", _current_chrom,":",i->first.first,"-",i->first.second, " has ",
              depth_i," read supporting. ","Intron at ", _current_chrom,":",j->first.first, "-",
              j->first.second, " has ", depth_j, " read supporting. ");
            }
         }
         if( depth_j / depth_i < kMinIsoformFrac){
            bad_intron_pos.push_back(j->first);
            if(_is_inspecting){
              LOG("Filtering overlapping intron by depth: ", _current_chrom,":",i->first.first,"-",i->first.second, " has ",
              depth_i," read supporting. ","Intron at ", _current_chrom,":",j->first.first, "-",
              j->first.second, " has ", depth_j, " read supporting. ");
            }
         }
       }
     }
   }
   sort(bad_intron_pos.begin(), bad_intron_pos.end());
   auto last = unique(bad_intron_pos.begin(), bad_intron_pos.end());
   bad_intron_pos.erase(last, bad_intron_pos.end());
   for(auto const& del: bad_intron_pos){
     auto to_del = intron_counter.find(del);
     assert(to_del != intron_counter.end());
     intron_counter.erase(to_del);
   }
//Filtering two: by small overhang supporting read proportion
//And at least two non small overhang supporting per intron

   for(auto i= intron_counter.cbegin(); i !=intron_counter.cend();){
     double total_read = i->second.total_junc_reads;
     double small_read = i->second.small_span_read;
//#ifdef DEBUG
//     cout<<small_read<<": ";
//     cout<<intron_counter[i].left<<" vs "<<intron_counter[i].right<<"\t"<<total_read<<endl;
//#endif
     if(total_read < kMinJuncSupport && !enforce_ref_models){
       if(_is_inspecting){
         LOG("Filtering intron at by overall read support: ", _current_chrom,":", i->first.first,"-",i->first.second,
            " has only ", total_read, " total read.");
       }
       i = intron_counter.erase(i);
       continue;
     }
     for(size_t k = i->first.first; k < i->first.second+1; ++k){
       intron_doc[k-cluster_left] += total_read;
     }
     if(small_read == total_read){
       i = intron_counter.erase(i);
       continue;
     }
     if(small_read < 1){
       ++i;
       continue;
     }
     double success = 2 * (double) kMaxSmallAnchor/_hit_factory->_reads_table._read_len_abs;;
     double prob_not_lt_observed = 1.0;
     double normal_mean = total_read * success;
     double normal_sd = sqrt(total_read * success*(1-success));
     double x = (small_read-0.5 - normal_mean)/normal_sd;
     prob_not_lt_observed = 1.0 - standard_normal_cdf(x);
     if(prob_not_lt_observed < kBinomialOverHangAlpha) {
       if(_is_inspecting){
         LOG("Filtering intron at by small anchor: ", _current_chrom,":", i->first.first,"-",i->first.second,
            " has ", small_read, " small overhang read vs ", total_read, " total read.");
       }
       i = intron_counter.erase(i);
       continue;
     }
     ++i;
   }

   //Filtering three: by comparing intron depth to exon depth
   for(auto i= intron_counter.begin(); i !=intron_counter.end();){
     uint start = i->first.first - cluster_left;
     uint end = i->first.second - cluster_left;
     assert(end > start);

     float avg_intron_doc = accumulate(intron_doc.begin()+start, intron_doc.begin()+end,0.0);
     avg_intron_doc /= (end-start);
     vector<float> exon_doc_dup(end-start+1);
     copy(exon_doc.begin()+start, exon_doc.begin()+end, exon_doc_dup.begin());

     //We were planning to calculate IRR which is a ratio of median depth of
     //read coverage on an intron to the number of read supporint this intorn
     //Right now the IRR hasn't been used yet.
     i->second.median_depth = getMedian(exon_doc_dup);


     float avg_intron_exonic_doc = accumulate(exon_doc_dup.begin(), exon_doc_dup.end(), 0.0);
     avg_intron_exonic_doc /= (end-start);
     if(avg_intron_exonic_doc != 0){
       if( avg_intron_doc / avg_intron_exonic_doc < kMinIsoformFrac){
         if(_is_inspecting){
            LOG("Filtering intron at by exonic coverage: ", _current_chrom, ":",i->first.first,"-",i->first.second,
            " averaged intron doc: ", avg_intron_doc, " vs averaged exonic doc on intron: ", avg_intron_exonic_doc, ".");
         }
         i = intron_counter.erase(i);
         continue;
       }
     }
     ++i;
     //cout << avg_intron_doc<<":"<<avg_intron_exonic_doc<<" left: "<<intron_counter[i].left \
         <<" right: "<<intron_counter[i].right<<endl;
   }
}


