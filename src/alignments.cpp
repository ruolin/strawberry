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

std::vector<vector<Contig>> Sample::assembleContig(const uint l, const uint r, const Strand_t strand, const vector<Contig>& hits) {
   std::vector<vector<Contig>> result;

   if (hits.empty()) {
      return result;
   }

   RefID ref_id = hits[0].ref_id();
   vector<float> exon_doc;
   IntronMap intron_counter;
   vector<GenomicFeature> exons;

   size_t s = r - l + 1;
   exon_doc.resize(s, 0);
   double avg_dep = 0.0;

//   for (const auto& h :hits) {
//     std::cerr<<h;
//   }
   avg_dep = compute_doc(l, r, hits, exon_doc, intron_counter, kMaxSmallAnchor);
   //cout<<avg_dep<<endl;

//#ifdef DEBUG
//   cout<<"cluster starts at: "<<l<<endl;
//   cout<<"cluster end at: "<<r<<endl;
//   cout<<"exon coverage:"<<endl;
//   for(auto i : exon_doc)
//      cout<<i;
//   cout<<endl;
//#endif

   if (avg_dep < kMinDepth4Locus) {
      return result;
   }

   filter_intron(this->_current_chrom, l, this->_hit_factory->_reads_table._read_len_abs, exon_doc, intron_counter);

   //local variables
   FlowNetwork flow_network;


   bool is_ok = flow_network.splicingGraph(ref_id, l, exon_doc, intron_counter, exons);
   if (!is_ok) return result;

#ifdef DEBUG
   std::cerr<<std::endl;
   for(auto e: exons){
      std::cerr<<"i exon: "<<e.left()<<"-"<<e.right()<<std::endl;
   }
   for(auto i:intron_counter){
       const IntronTable &intron = i.second;
       std::cerr<<"i intron: "<<intron.left<< "-" << intron.right<<std::endl;
   }
#endif
   vector<vector<GenomicFeature>::const_iterator> bundles;
   bundles.push_back(exons.cbegin());
   for (auto it = exons.cbegin(); it + 1 < exons.cend(); ++it) {
      auto ij = it + 1;
      if (it->right() + 1 == ij->left()) continue;

      bool no_connect = true;
      //search forward
      for (auto search = it + 1; search < exons.cend(); ++search) {
         if (intron_counter.find(make_pair<uint,uint>(it->right() + 1, search->left() - 1)) != intron_counter.end()) {
            no_connect = false;
         }
      }
      //search reverse
      for (auto search = it; search >= exons.cbegin(); --search) {
         if (intron_counter.find(make_pair<uint,uint>(search->right() + 1, ij->left() - 1)) != intron_counter.end()) {
            no_connect = false;
         }
      }
      if (no_connect) {
         LOG(INFO)<<ref_id<<":"<<l<<"-"<<r<<"refine graph by spliting";
         bundles.push_back(ij);
         bundles.push_back(ij);
      }
   }
   bundles.push_back(exons.cend());
   for (size_t i = 0; i < bundles.size(); i += 2) {
      vector<GenomicFeature> refined_exons(bundles[i], bundles[i+1]);
      std::map<std::pair<uint,uint>, IntronTable> refined_introns;
      for (auto const& intron : intron_counter) {
         uint right_bound = bundles[i+1] == exons.cend() ? std::numeric_limits<uint>::max() : bundles[i+1]->left();
         if (intron.first.second < right_bound && intron.first.first > bundles[i]->right()) {
            refined_introns.emplace(intron.first, intron.second);
         }
      }
      auto txs = runFlowAlgorithm(strand, hits, refined_introns, refined_exons);
      result.push_back(vector<Contig>(txs.begin(), txs.end()));
   }
   return result;
}

vector<Contig> Sample::runFlowAlgorithm(const Strand_t& strand, const vector<Contig>& hits,
                                      const std::map<std::pair<uint,uint>, IntronTable> &intron_counter,
                                      const std::vector<GenomicFeature> &exons) {
   FlowNetwork flow_network;
   Graph::NodeMap<const GenomicFeature *> node_map(flow_network._g);
   Graph::ArcMap<int> cost_map(flow_network._g);
   Graph::ArcMap<int> min_flow_map(flow_network._g);
   vector<vector<Graph::Arc>> path_cstrs;
   vector<vector<GenomicFeature>> assembled_feats;
   vector<vector<size_t>> constraints;


   vector<Contig> result;
   bool stat = flow_network.createNetwork(hits, exons, intron_counter, node_map, cost_map, min_flow_map, path_cstrs);
   if (!stat) {
      //std::cerr<<"warning: cannot create network!\n";
      return result;
   }

   bool stat2 = flow_network.solveNetwork(node_map, exons, path_cstrs, cost_map, min_flow_map, assembled_feats);
   if (!stat2) {
      std::cerr<<"warning: cannot slove network!\n";
      return result;
   }
   RefID ref_id = hits[0].ref_id();
   return assemble_2_contigs(assembled_feats, ref_id, strand);
}
/*
 * Global utility functions end:
 */

HitCluster::HitCluster():
   _leftmost(UINT_MAX),
   _rightmost(0),
   //_plus_strand_num_hits(0),
   //_minus_strand_num_hits(0),
   //_first_encounter_strand(Strand_t::StrandUnknown),
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

void HitCluster::addRefContig(const Contig& contig)
{
   if(_ref_id != -1){
     assert(_ref_id == contig.ref_id());
   }
   else{
     _ref_id = contig.ref_id();
   }
   //if (gene_id().empty()) gene_id() = contig.parent_id();
   //else
   if (gene_id() != contig.parent_id()) {
      //std::cerr<<gene_id()<<" and "<<contig.parent_id()<<" overlaps!\n";
      return;
   }
   _leftmost = min(_leftmost, contig.left());
   _rightmost = max(_rightmost, contig.right());
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
   return _hits.size();
}

int HitCluster::len() const{
   return _rightmost - _leftmost + 1;
}

Strand_t HitCluster::ref_strand() const{
   assert(!_ref_mRNAs.empty());
   Strand_t strand = _ref_mRNAs[0].strand();
   for(auto &i: _ref_mRNAs){
     assert(i.strand() == strand);
   }
   return strand;
}

Strand_t HitCluster::guessStrand() const{
   int max_count = INT_MIN;
   Strand_t best_strand = Strand_t::StrandUnknown;
   if (_strand_intron.find(Strand_t::StrandPlus) != _strand_intron.end()) {
      for (auto const& i : _strand_intron.at(Strand_t::StrandPlus)) {
         //std::cerr<<"+, "<<i.first<<":"<<i.second<<std::endl;
         if (i.second > max_count) {
            max_count = i.second;
            best_strand = Strand_t::StrandPlus;
         }
      }
   }
   if (_strand_intron.find(Strand_t::StrandMinus) != _strand_intron.end()) {
      for (auto const& i : _strand_intron.at(Strand_t::StrandMinus)) {
         //std::cerr<<"-, "<<i.first<<":"<<i.second<<std::endl;
         if (i.second > max_count) {
            max_count = i.second;
            best_strand = Strand_t::StrandMinus;
         }
      }
   }
   return best_strand;
}

//void HitCluster::reset(uint old_left,
//                  uint old_right,
//                  RefID old_ref_id,
//                  int old_plus_hit,
//                  int old_minus_hit,
//                  Strand_t old_strand,
//                  double old_mass){
//   _leftmost = old_left;
//   _rightmost = old_right;
//   _ref_id = old_ref_id;
//   _plus_strand_num_hits = old_plus_hit;
//   _minus_strand_num_hits = old_minus_hit;
//   _first_encounter_strand = old_strand;
//   _raw_mass = old_mass;
//}

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
   if (hit._left_read && hit._left_read->contains_splice()) {
      vector<GenomicFeature> gfs;
      if (readhit_2_genomicFeats(hit.left_read_obj(), gfs)) {
         for (auto const & gf : gfs) {
            if (gf._match_op._code == Match_t::S_INTRON) {
               _strand_intron[hit._left_read->strand()][gf]++ ;
            }
         }
      }
   }
   if (hit._right_read && hit._right_read->contains_splice()) {
      vector<GenomicFeature> gfs;
      if (readhit_2_genomicFeats(hit.right_read_obj(), gfs)) {
         for (auto const & gf : gfs) {
            if (gf._match_op._code == Match_t::S_INTRON) {
               _strand_intron[hit._right_read->strand()][gf]++;
            }
         }
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
   uint orig_left = _leftmost;
   uint orig_right = _rightmost;
   RefID orig_ref_id = _ref_id;
   uint hit_left = hit->left();
   uint hit_right = hit->right();
   Strand_t hit_strand = hit->strand();
   RefID hit_ref_id = hit->ref_id();
   uint hit_partner_pos = hit->partner_pos();
   ReadID hit_id = hit->read_id();

   if(extend_by_hit){
     _leftmost = min(_leftmost, hit_left);
     _rightmost = max(_rightmost, hit_right);
   }
   if(extend_by_partner && hit_partner_pos != 0 && hit->partner_ref_id() == _ref_id){
     if((int)hit_partner_pos - (int)hit_left < kMaxIntronLength){
       _rightmost = max(max(_rightmost, hit->right()), hit->partner_pos());
     }
   }

   // Double check. This is only useful when called in Sample::mergeCluster()
//   if(hit_lt_cluster(*hit, *this, kMaxOlapDist)){
//     reset(orig_left, orig_right, orig_ref_id);
//     return false;
//   }

   if(abs((int)hit_right - (int)hit_left) > kMaxFragSpan){
     reset(orig_left, orig_right, orig_ref_id);
      if (!NO_LOGGING) {
         LOG(WARNING)<<"Hit start at "<<hit_left<< "  is longer than max gene length<< skipping";
      }
     return false;
   }

   if(_ref_id == -1){
     if(hit_ref_id != -1)
       _ref_id = hit_ref_id;
   } else{
     assert(_ref_id == hit_ref_id);
   }

   if(hit->is_singleton() || hit->partner_ref_id() != _ref_id){
   //if(hit->is_singleton()){
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


       if(hit->partner_pos() > hit->left()){
         if(hit->reverseCompl())
         {
            if (!NO_LOGGING) {
               LOG(WARNING)<<"Possible wrong read orientation at chr: "<< hit->ref_id()<< " for read start at "<< hit->left()<< " and his partner at "<<hit->partner_pos();
            }
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
            if (!NO_LOGGING) {
               LOG(WARNING)<<"Possible wrong read orientation at chr: "<< hit->ref_id()<< " for read start at "<< hit->left()<< " and his partner at "<<hit->partner_pos();
            }
         }

         PairedHit open_hit(nullptr, hit);
         unordered_map<ReadID, list<PairedHit>>::iterator ins_pos;
         list<PairedHit> chain;
         chain.push_back(move(open_hit));
         bool status;
         tie(ins_pos, status) = _open_mates.insert(make_pair(hit_id, chain));
         assert(status);
       }
       else{ // hit and its partner start at the some position
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
       if(hit->partner_pos() > hit->left()){
         PairedHit open_hit(hit, nullptr);
         iter_open->second.push_back(open_hit);
       }
       else if(hit->partner_pos() < hit->left()){
         PairedHit open_hit(nullptr, hit);
         iter_open->second.push_back(open_hit);
       }
       else{ //read and its partner starts at some position
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
     assert(false);
   }
   if(_hits.empty())
     return 0;

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

void HitCluster::setBoundaries(){
/*
 * Call this after collapseHits
 */
   if(enforce_ref_models && hasRefmRNAs()){
     _leftmost = INT_MAX;
     _rightmost = 0;
     for(auto r : _ref_mRNAs){
       _leftmost = min(_leftmost, r.left());
       _rightmost = max(_rightmost, r.right());
     }
   }
}


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
   int plus_count = 0;
   int minus_count = 0;
   if (_strand_intron.find(Strand_t::StrandPlus) != _strand_intron.end()) {
      for (auto const& i : _strand_intron.at(Strand_t::StrandPlus)) {
         plus_count += i.second;
      }
   }
   if (_strand_intron.find(Strand_t::StrandMinus) != _strand_intron.end()) {
      for (auto const& i : _strand_intron.at(Strand_t::StrandMinus)) {
         minus_count += i.second;
      }
   }

   int minor = std::min(plus_count, minus_count);
   int major = std::max(plus_count, minus_count);
#ifdef DEBUG
   std::cerr<<"minor iso junc: "<<minor<<std::endl;
   std::cerr<<"major iso junc: "<<major<<std::endl;
#endif
   int scale = 1;
   if (minor > (int) scale * major * kMinIsoformFrac) return true;
   return false;
}

bool Sample::load_chrom_fasta(RefID seq_id)
{
   /*
   * Loading reference FASTA sequence without gene model
   */
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
   cerr<<"Has loaded transcripts from "<<gseqs.size()<<" Chromosomes/Scaffolds"<<endl;
   /*
   * Parse transcripts in GffTree to a vector of Contig objects.
   */

   //sort gseqs accroding to the observation order in ref_table
   // or if ref_table is empty, initialize it according to gseqs
   //ref_id in ref_table start from 0.
   if(rt.size() == 0){
     for(uint i=0; i<gseqs.size(); ++i){
       rt.set_id(gseqs[i]->_g_seq_name);
     }
   } else{
     for(uint i = 0; i<gseqs.size(); ++i){
       int idx = gseqs[i]->get_gseq_id();
       int ref_table_id = rt.get_id(gseqs[i]->_g_seq_name);
       if (ref_table_id == -1) {
           cerr<<"Warning: the gff file contains seq name "<<gseqs[i]->_g_seq_name<<" which is not found in the bam file"<<endl;
       }
       else if(idx != ref_table_id ){
         cerr<<"Warning: Sam file and Gff file are not sorted in the same order!\n";
         cerr<<"set gff chrom id to "<<ref_table_id<<endl;
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

//     if(_prev_hit_ref_id != -1){
//       if(_prev_hit_ref_id == readin.ref_id() && _prev_hit_pos > readin.left())
//       {
//         const string cur_chr_name = _hit_factory->_ref_table.ref_real_name(readin.ref_id());
//         const string last_chr_name = _hit_factory->_ref_table.ref_real_name(_prev_hit_ref_id);
//         cerr<<"BAM file not sort correctly!\n";
//         cerr<<"The current position is: "<<cur_chr_name<<":"<<readin.left();
//         cerr<<" and previous position is: "<<last_chr_name<<":"<<_prev_hit_pos;
//         cerr<<endl;
//       }
//     }
//
//     _prev_hit_ref_id = readin.ref_id();
//     _prev_hit_pos = readin.left();
     break;
   }
   _hit_factory->_reads_table._read_len_abs = max(_hit_factory->_reads_table._read_len_abs, readin.read_len());
   return raw_mass;
}

void Sample::rewindHit()
{
   //double mass = rh.mass();
   _hit_factory->undo_hit();
   //return mass;
}

void Sample::rewindHits(uint64_t pos)
{
   //double mass = rh.mass();
   _hit_factory->return2Pos(pos);
   //return mass;
}

int Sample::addRef2Cluster(HitCluster &cluster_out){
   if(_refmRNA_offset >=  _ref_mRNAs.size()) {
      //std::cerr<<" all transcripts processed!\n";
     _has_load_all_refs = true;
     return 0;
   }

   // add first rna
   //cout<<"offset: "<<_refmRNA_offset<<endl;
   //cout<<"_ref_mRNAs size: "<<_ref_mRNAs[0].parent_id()<<endl;

   cluster_out.gene_id() = _ref_mRNAs[_refmRNA_offset].parent_id();
   //std::cerr<<"a: "<<cluster_out.gene_id()<<std::endl;
   cluster_out.addRefContig(_ref_mRNAs[_refmRNA_offset++]);
   if(_refmRNA_offset >= _ref_mRNAs.size()){
     _has_load_all_refs = true;
     return 1;
   }

   if (!cluster_out.gene_id().empty()) {
      while (_refmRNA_offset < _ref_mRNAs.size() && _ref_mRNAs[_refmRNA_offset].parent_id() == cluster_out.gene_id()) {
         cluster_out.addRefContig(_ref_mRNAs[_refmRNA_offset++]);
      }
      size_t mark_next_gene = _refmRNA_offset;
      //continue search a few forward
      int over = 0;
      while (++_refmRNA_offset < _ref_mRNAs.size() && over++ < 10) {
         if (_ref_mRNAs[_refmRNA_offset].parent_id() == cluster_out.gene_id()) {
            cluster_out.addRefContig(_ref_mRNAs[_refmRNA_offset]);
         }
      }
      _refmRNA_offset = mark_next_gene;
   }
   else {
      // add the rest if overlapped with first
      size_t i = 0;
      while (i < cluster_out._ref_mRNAs.size()) {
         const Contig &ref = cluster_out._ref_mRNAs[i];
         if (Contig::overlaps_directional(ref, _ref_mRNAs[_refmRNA_offset])) {
            cluster_out.addRefContig(_ref_mRNAs[_refmRNA_offset++]);
            if (_refmRNA_offset >= _ref_mRNAs.size()) {
               _has_load_all_refs = true;
               return cluster_out._ref_mRNAs.size();
            }
            i = 0;
         } else {
            ++i;
         }
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

     if(!_hit_factory->recordsRemain()){
       return clusterOut.size();
     }

     if(new_hit->ref_id() > next_ref_start_ref ||
      (new_hit->ref_id() == next_ref_start_ref && new_hit->right() >= next_ref_start_pos)){
       rewindHit();
       return clusterOut.size();
     }

     if(clusterOut.ref_id() == -1){ // add first hit

       clusterOut.addOpenHit(new_hit, true, true);
       clusterOut.addRawMass(mass);
     } else { //add the rest
       if(hit_lt_cluster(*new_hit, clusterOut, kMaxOlapDist)){
         // should never reach here
         std::cerr<<"It appears that SAM/BAM not sorted!\n";
         continue;
       }
       if(hit_gt_cluster(*new_hit, clusterOut, kMaxOlapDist)){
         // read has gone to far.
         rewindHit();
         break;
       }
       clusterOut.addOpenHit(new_hit, true, true);
       clusterOut.addRawMass(mass);
     }
   }
   return clusterOut.size();
}

int Sample::nextClusterRefDemand(HitCluster &clusterOut){
   // if assembly step was run. hasLoadRefmRNAs() will return true.

   if (!hasLoadRefmRNAs()) {
     std::cerr<<"if you use --no-assembly option, you must provide gff file through -g option!"<<std::endl;
     assert(false);
   }
   if (!_hit_factory->recordsRemain()) {
      return -1;
   }
   int num_added_refmRNA = addRef2Cluster(clusterOut);
   //cout<<num_added_refmRNA<<" added mRNA"<<endl;
   if (num_added_refmRNA == 0) {
     return -1;
   }
   int64_t first_pos = _hit_factory->getCurrPos();
//   int counter = 0;
   while (true) {
     if(!_hit_factory->recordsRemain()){
       break;
     }
     ReadHitPtr new_hit(new ReadHit());
     double mass = next_valid_alignment(*new_hit);
     if (hit_lt_cluster(*new_hit, clusterOut, 0)) {  //hit hasn't read this region

     } else if (hit_gt_cluster(*new_hit, clusterOut, 0)) {
       rewindHit();
       break;
     } else if (new_hit->strand() != Strand_t ::StrandUnknown && new_hit->strand() != clusterOut.ref_strand()) {
     }
     else {
       clusterOut.addOpenHit(new_hit, false, false);
       clusterOut.addRawMass(mass);
     }
   }  //end while loop

   //roger
//   if (clusterOut.size() > 0) {
//      std::cout<<num_added_refmRNA<<" transcript added\n";
//   }
   return clusterOut.size();
}

void Sample::preProcess(FILE *log) {
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

       if(hit_lt_cluster(*new_hit, clusterOut, kMaxOlapDist)){ // hit hasn't reach reference region
         rewindHit();
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

       if(hit_gt_cluster(*new_hit, clusterOut, kMaxOlapDist)){ // read has gone too far.
         rewindHit();
         break;
       }

       clusterOut.addOpenHit(new_hit, false, false);
       clusterOut.addRawMass(mass);
       if (!_hit_factory->recordsRemain()) return clusterOut.size();
     } // end while loop
   } // end loadRefmRNAs
   return clusterOut.size();
}

void Sample::mergeClusters(HitCluster & last, HitCluster &cur){
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
   cluster->setBoundaries(); // set boundaries if reference exist.
}

void Sample::fragLenDist(const RefSeqTable &ref_t,
               const std::vector<Contig> &transcripts,
               const shared_ptr<HitCluster> cluster,
               FILE *plogfile) {

   if (transcripts.empty()) {
      LOG(WARNING) << "no reference transcripts are found";
      return;
   }
   _total_mapped_reads += (int) cluster->weighted_mass();

   vector<Contig> hits;
   for (auto r = cluster->_uniq_hits.cbegin(); r != cluster->_uniq_hits.cend(); ++r) {
     Contig hit(*r);
     if (hit.ref_id() != -1) {
       hits.push_back(hit);
     }
   }

   //if (transcripts.size() == 1) {
   for (size_t h = 0; h < hits.size(); ++h) {
     int counter = 0;
     size_t mark = 0;
     for (size_t t = 0; t < transcripts.size(); ++t) {
        //if (hits[h].is_single_read()) continue;
        if (Contig::is_compatible(hits[h], transcripts[t])) {
           ++counter;
           mark = t;
        }
     } //end for
     if (counter == 1) {
        double frag_len = Contig::exonic_overlaps_len(transcripts[mark],
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
     }
#else
         _hit_factory->_reads_table._frag_dist.push_back(frag_len);
#endif
   }// end for

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
   vector<Contig> hits;

   if (cluster->size() < kMinReadForAssemb) {
//#if ENABLE_THREADS
//      if (use_threads) decr_pool_count();
//#endif
      return assembled_transcripts;
   }

   cluster->_strand = cluster->guessStrand();
   Strand_t second_strand = Strand_t::StrandUnknown;
   bool see_both_strand = cluster->see_both_strands();

   if(see_both_strand) {
      switch(cluster->_strand) {
         case Strand_t::StrandPlus:
            second_strand = Strand_t::StrandMinus;
            break;
         case Strand_t::StrandMinus:
            second_strand = Strand_t::StrandPlus;
            break;
      }
   }

   vector<Contig> hits_on_other_strands;
   for (auto r = cluster->_uniq_hits.cbegin(); r != cluster->_uniq_hits.cend(); ++r) {
      Contig hit(*r);
      if (hit.ref_id() != -1) {
         if (hit.strand() == cluster->_strand || hit.strand() == Strand_t::StrandUnknown) {
            hits.push_back(hit);
         }
         if (see_both_strand) {
            if (hit.strand() == second_strand || hit.strand() == Strand_t::StrandUnknown) {
               hits_on_other_strands.push_back(hit);
            }
         }
      }
   }

   if (cluster->hasRefmRNAs() && utilize_ref_models ) {
      for (const auto& i: cluster->_ref_mRNAs) {
         if (i.strand() == cluster->_strand) {
            hits.push_back(i);
            hits.back()._is_ref = true;
         } else {
            hits_on_other_strands.push_back(i);
            hits_on_other_strands.back()._is_ref = true;
         }
      }
   }
   sort(hits.begin(), hits.end());
   auto result = this->assembleContig(cluster->left(), cluster->right(), cluster->strand(), hits);
   for (auto& transcripts: result) {
      AddTranscripts(transcripts, assembled_transcripts);
   }

   if (see_both_strand) {
      sort(hits_on_other_strands.begin(), hits_on_other_strands.end());
      auto result = this->assembleContig(cluster->left(), cluster->right(), second_strand, hits_on_other_strands);
      for (auto& transcripts: result) {
         AddTranscripts(transcripts, assembled_transcripts);
      }
   }

   this->fragLenDist(ref_t, assembled_transcripts, cluster, plogfile);

   return assembled_transcripts;
}

void Sample::AddTranscripts(std::vector<Contig>& transcripts, std::vector<Contig>& result) {
   ++_num_cluster;
   int tid=0;
   for (Contig& asmb: transcripts) {
      ++tid;
      asmb.parent_id() = "gene."+to_string(_num_cluster);
      asmb.annotated_trans_id("transcript." + to_string(_num_cluster) + "." + to_string(tid));
   }
   result.insert(result.end(), transcripts.begin(), transcripts.end());
}

void Sample::quantifyCluster(const RefSeqTable &ref_t, const shared_ptr<HitCluster> cluster,
                 const vector<Contig> &assembled_transcripts, FILE *pfile, FILE *plogfile, FILE *fragfile) const {


   LocusContext est(*this, plogfile, cluster, assembled_transcripts);

   //est.set_empirical_bin_weight(iso_2_bins_map, iso_2_len_map, cluster->collapse_mass(), exon_bin_map);
   //est.calculate_raw_iso_counts(iso_2_bins_map, exon_bin_map);
   bool success = est.estimate_abundances();

   if(success){
      //cout<<"assembled transcripts size: "<<assembled_transcripts.size()<<endl;
#if ENABLE_THREADS
     if(use_threads) {
        out_file_lock.lock();
     }
#endif
     cerr << ref_t.ref_real_name(cluster->ref_id()) << "\t" << cluster->left() << "\t" << cluster->right()
          << " finishes abundances estimation" << endl;

     for (const auto &iso: est.transcripts()) {
        iso._contig.print2gtf(pfile, _hit_factory->_ref_table, iso._FPKM_s,
                                 iso._frac_s, iso._gene_str, iso._isoform_str);
     }

      if (fragfile != NULL) {
         //if (est.num_transcripts() > 1 && est.num_exon_bins() > 1) {
            printContext(est, cluster, _fasta_getter, fragfile);
         //}
      }
   }
#if ENABLE_THREADS
   if(use_threads) {
     out_file_lock.unlock();
     decr_pool_count();
   }
#endif
}


void Sample::printContext(const LocusContext& est, const shared_ptr<HitCluster> cluster,
                          const std::shared_ptr<FaSeqGetter> & fa_getter, FILE *fragfile) const {
   /* Print locus coordinates*/
   map<set<pair<uint,uint>>, uint> eb_count_map;
   map<set<pair<uint,uint>>, vector<double>> eb_prob_map;
   for (auto r = cluster->uniq_hits().cbegin(); r != cluster->uniq_hits().cend(); ++r) {
      Contig hit = Contig(*r);
      if (hit.ref_id() == -1) continue;
      auto eb = est.get_frag_info(hit);
      if (!eb.first.empty()) {
         ++eb_count_map[eb.first];
         eb_prob_map[eb.first] = eb.second;
      }
   }

   uint sum = 0;
   for (const auto& eb: eb_count_map) {
      sum += eb.second;
   }


   for (auto it = eb_prob_map.cbegin(); it != eb_prob_map.cend(); ++it) {
      vector<string> info;
      info.push_back(est.sample_name());
      info.push_back(to_string(total_mapped_reads()));

      info.push_back(est.gene_name());
      info.push_back(to_string(sum));

      string iso_names;
      for (const auto& tn : est.transcript_names()) {
         iso_names += tn;
         iso_names += ",";
      }
      iso_names.pop_back();
      info.push_back(iso_names);

      string cond_prop;
      for (const double prob:it->second) {
         cond_prop += to_string_with_precision(prob, 12);
         cond_prop += ",";
      }
      cond_prop.pop_back();
      info.push_back(cond_prop);

      string class_prop;
      for (auto iso = est.transcripts().cbegin(); iso != est.transcripts().cend(); ++iso) {
         class_prop += iso->_frac_s;
         class_prop += ",";
      }
      class_prop.pop_back();
      info.push_back(class_prop);

      string coords;
      for (auto const& c: it->first) {
         coords += "[";
         coords += to_string(c.first);
         coords += "-";
         coords += to_string(c.second);
         coords += "]";
      }
      info.push_back(coords);
      info.push_back(to_string(eb_count_map[it->first]));

      if (BIAS_CORRECTION) {
         auto seq = ExonBin::bin_dnaseq(it->first, fa_getter);
         auto gcratio = Kmer<string>::GCRatio(seq.begin(), seq.end());
         auto entropy = Kmer<string>::Entropy(seq, 6); // hexmer entropy
         auto highgc2080 = Kmer<string>::HighGCStrech(seq.begin(), seq.end(), 20, 0.8);
         auto highgc2090 = Kmer<string>::HighGCStrech(seq.begin(), seq.end(), 20, 0.9);
         auto highgc4080 = Kmer<string>::HighGCStrech(seq.begin(), seq.end(), 40, 0.8);
         auto highgc4090 = Kmer<string>::HighGCStrech(seq.begin(), seq.end(), 40, 0.9);
         info.push_back(to_string(gcratio));
         info.push_back(to_string(entropy));
         info.push_back(to_string(highgc2080));
         info.push_back(to_string(highgc2090));
         info.push_back(to_string(highgc4080));
         info.push_back(to_string(highgc4090));
      }
      if (fragfile != NULL) pretty_print(fragfile, info, "\t");
   }
}

void Sample::addAssembly(const std::vector<Contig>& assembs) {

#if ENABLE_THREADS
   if (use_threads) thread_pool_lock.lock();
#endif

   if (!assembs.empty()) std::move(assembs.begin(), assembs.end(), std::back_inserter(_assembly));

#if ENABLE_THREADS
   if(use_threads) thread_pool_lock.unlock();
#endif
}

void Sample::assembleSample(FILE *plogfile)
/*
 *  First run-through to calculate fragment distribution FD
 * */
{
   const RefSeqTable & ref_t = _hit_factory->_ref_table;
   _num_cluster = 0;

   RefID current_ref_id = std::numeric_limits<int>::max();

   while(true){
     shared_ptr<HitCluster> cur_cluster (new HitCluster());
     if(-1 == nextCluster_refGuide(*cur_cluster)){
       break;
     }
     if(cur_cluster->ref_id() == -1){
       continue;
     }

//Begin loading ref seqs
     if(current_ref_id != cur_cluster->ref_id()){
       current_ref_id = cur_cluster->ref_id();
     }
     if(_current_chrom != ref_t.ref_real_name(cur_cluster->ref_id())){
       _current_chrom = ref_t.ref_real_name(cur_cluster->ref_id());
     }
//End loading ref seqs
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
            finalizeCluster(cur_cluster, true);
            vector<Contig> asmb = this-> assembleCluster(ref_t, cur_cluster, plogfile);
            this->addAssembly(asmb);
            --curr_thread_num;
            });
       //thread worker{&Sample::finalizeAndassembleCluster, this, ref_t, cur_cluster, NULL, NULL};
       worker.detach();
     }else{
       finalizeCluster(cur_cluster, true);
       vector<Contig> asmb = assembleCluster(ref_t, cur_cluster, plogfile);
       this->addAssembly(asmb);
     }
#else
   finalizeCluster(cur_cluster, true);
   assembleCluster(ref_t, cur_cluster, plogfile);
   this->addAssembly(asmb);
#endif
     // _total_mapped_reads += (int) cur_cluster->raw_mass();
     //cout<<"weighted cluster mass: "<<cur_cluster->_weighted_mass<<endl;
     fprintf(plogfile, "Inspect gene: %s:%d-%d\n", ref_t.ref_real_name(cur_cluster->ref_id()).c_str(), cur_cluster->left(), cur_cluster->right());
     fprintf(plogfile, "Has inspected %d reads\n", (int)_total_mapped_reads);
     cur_cluster = move(cur_cluster);
   }

//make sure all threads have finished
#if ENABLE_THREADS
   if(use_threads){
     while(true){
       if(curr_thread_num==0){
         break;
       }
       this_thread::sleep_for(chrono::milliseconds(3));
     }
   }
#endif
}

int Sample::total_mapped_reads() const
{
   return _total_mapped_reads;
}

void Sample::procSample(FILE *pfile, FILE *plogfile, FILE *fragfile)
{
/*
 * The major function which calls nextCluster() and finalizes cluster and
 * assemble each cluster-> Right now only nextCluster_refGuide() is implemented.
 * if no reference mRNA than nextCluster_refGuide will call nextCluster_denovo()
 */
   NO_LOGGING = true;
   _hit_factory->reset();
   reset_refmRNAs();
   const RefSeqTable & ref_t = _hit_factory->_ref_table;
   int current_ref_id = INT_MAX;

   if (fragfile != NULL) {
      //if (est.num_transcripts() > 1 && est.num_exon_bins() > 1) {
      std::vector<string> header = {"sample", "sample_frag_count", "gene_id", "gene_frag_count",
      "transcripts", "conditonal_probabilites", "class_probabilites", "path_symbol", "path_count",
      "path_gc_content", "path_hexmer_entropy", "gc_strecth_0.8_20", "gc_strecth_0.9_20", "gc_strecth_0.8_40",
      "gc_strecth_0.9_40"};
      if (fragfile != NULL) pretty_print(fragfile, header, "\t");
      //}
   }

   while(true){
     //++_num_cluster;
     shared_ptr<HitCluster> cluster (new HitCluster());
     if(-1 == nextClusterRefDemand(*cluster)){
       break;
     }
     if(cluster->ref_id() == -1) continue;

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
            finalizeCluster(cluster, true);
            this->quantifyCluster(ref_t, cluster, cluster->ref_mRNAs(), pfile, plogfile, fragfile);
       });
       worker.detach();
     }else {
       finalizeCluster(cluster, true);
       this->quantifyCluster(ref_t, cluster, cluster->ref_mRNAs(), pfile, plogfile, fragfile);
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


double compute_doc(const uint left, const uint right,
                   const vector<Contig> & hits,
                   vector<float> &exon_doc,
                   IntronMap &intron_counter, uint smallOverHang)
{

   assert(right > left);
   for(size_t i = 0; i<hits.size(); ++i){
     const vector<GenomicFeature> & g_feats = hits[i]._genomic_feats;
     for(size_t j = 0; j<g_feats.size(); ++j){
       const GenomicFeature & gf = g_feats[j];
       if( gf._match_op._code == Match_t::S_MATCH){
         size_t l  = std::max(left, gf.left());
         size_t r = std::min(gf.right(), right);
         for(size_t p = l; p < r+1; ++p){
            exon_doc[p-left] += hits[i].mass();
         }
       }
       else if( gf._match_op._code == Match_t::S_INTRON){
         if(gf.left() < left || gf.right() > right)
            continue;
         if (gf.left() >= gf.right()) continue;
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

void filter_intron(const std::string& current_chrom, const uint cluster_left,
                   const uint read_abs_len, const vector<float> &exon_doc, IntronMap& intron_counter)
{
   // bad_intron_pos for indexes for intron to be drop off in intron_counter
   vector<float> intron_doc(exon_doc.size(),0.0);

//Filtering one: by overlapping with better intron
   vector<pair<uint, uint>> bad_intron_pos;
   for(auto i = intron_counter.cbegin(); i != intron_counter.cend(); ++i){
     for(auto j = next(i); j != intron_counter.cend(); ++j){
        int scale = -1;
        if(IntronTable::overlap(i->second, j->second)){
           scale = 1;
           if (!IntronTable::contains_or_is_contained(i->second, j->second)) {
             scale = 10;
           }
        }

        float depth_i = i->second.total_junc_reads;
        float depth_j = j->second.total_junc_reads;
        std::pair<uint, uint> bad;
        float min_junc = 0;
        if (depth_j < depth_i) {
           min_junc = depth_j;
           bad = j->first;
        } else {
           min_junc = depth_i;
           bad = i->first;
        }
        if( min_junc / (depth_i + depth_j) < kMinIsoformFrac * scale){
           bad_intron_pos.push_back(bad);
           if (!NO_LOGGING) {
              LOG(INFO)<<"Filtering overlapping intron by depth: "<< current_chrom<<":"<<i->first.first<<"-"<<i->first.second<< " has "<<
                       depth_i<<" read supporting. "<<"Intron at "<< current_chrom<<":"<<j->first.first<< "-"<<
                       j->first.second<< " has "<< depth_j<< " read supporting. ";
           }
        }
     } // end for
   } //end for
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
     if (total_read < kMinJuncSupport && !enforce_ref_models) {
        if (!NO_LOGGING) {
           LOG(INFO)<<"Filtering intron at by overall read support: "<< current_chrom<<":"<< i->first.first<<"-"<<i->first.second<<
                    " has only "<< total_read<< " total read.";
        }
        i = intron_counter.erase(i);
        continue;
     }
     if(i->first.second - i->first.first > LongJuncLength && total_read < kMinSupportForLongJunc && !enforce_ref_models){
        if (!NO_LOGGING) {
           LOG(INFO)<<"Filtering long intron at by overall read support: "<< current_chrom<<":"<< i->first.first<<"-"<<i->first.second<<
                    " has only "<< total_read<< " total read.";
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
     double success = 2 * (double) kMaxSmallAnchor / read_abs_len;
     double prob_not_lt_observed = 1.0;
     double normal_mean = total_read * success;
     double normal_sd = sqrt(total_read * success*(1-success));
     double x = (small_read-0.5 - normal_mean)/normal_sd;
     prob_not_lt_observed = 1.0 - standard_normal_cdf(x);
     if(prob_not_lt_observed < kBinomialOverHangAlpha) {
        if (!NO_LOGGING) {
           LOG(INFO)<<"Filtering intron at by small anchor: "<< current_chrom<<":"<< i->first.first<<"-"<<i->first.second<<
                    " has "<< small_read<< " small overhang read vs "<< total_read<< " total read.";
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
         if (!NO_LOGGING) {
            LOG(INFO)<<"Filtering intron at by exonic coverage: "<< current_chrom<< ":"<<i->first.first<<"-"<<i->first.second<<
                      " averaged intron doc: "<< avg_intron_doc<< " vs averaged exonic doc on intron: "<< avg_intron_exonic_doc<< ".";
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



