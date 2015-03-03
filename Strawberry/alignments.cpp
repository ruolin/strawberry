/*
 * alignment.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: ruolin
 */

#include "alignments.h"
#include "fasta.h"
#include <algorithm>
#include <iterator>
#include <random>
#ifdef DEBUG
   #include <iostream>
#endif
using namespace std;

/*
 * Global utility function:
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

bool hit_complete_within_cluster(const PairedHit& hit, const HitCluster& cluster, uint olap_radius){
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
   return _hits.size();
}

Strand_t HitCluster::ref_strand() const{
   assert(!_ref_mRNAs.empty());
   Strand_t strand = _ref_mRNAs[0]->strand();
   for(auto &i: _ref_mRNAs){
      assert(i->strand() == strand);
   }
   return strand;
}

bool HitCluster::addHit(const PairedHit &hit){
   if(_final){
      return false;
   }
   assert(_ref_id == hit.ref_id());
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


bool HitCluster::addOpenHit(ReadHitPtr hit, bool extend_by_hit, bool extend_by_partner)
{
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
   if(extend_by_partner){
      _rightmost = max(max(_rightmost, hit->right()), hit->partner_pos());
   }

   // Double check. This is only useful when called in ClusterFactory::mergeCluster()
   if(hit_lt_cluster(*hit, *this, ClusterFactory::_kMaxOlapDist)){
      return false;
   }

   if(abs((int)hit_right - (int)hit_left) > _kMaxGeneLen){

      LOG_WARN("Hit start at ",hit_left, "  is longer than max gene length, skipping");
      return false;
   }
   _ref_id = hit_ref_id;
   if(hit->is_singleton()){
      PairedHit ph(hit, nullptr);
      addHit(ph);
   }
   else{
      unordered_map<ReadID, PairedHit>::iterator iter_open = _open_mates.find(hit_id);
      if( iter_open == _open_mates.end()){
         PairedHit open_hit(hit, nullptr);
         unordered_map<ReadID, PairedHit>::iterator ins_pos;
         bool status;
         tie(ins_pos, status) = _open_mates.insert(make_pair(hit_id, move(open_hit)));
         assert(status);
      } else{

         if(iter_open->second.ref_id() != hit_ref_id)
            return false;

         bool strand_agree = iter_open->second.strand() == hit_strand ||
         iter_open->second.strand() == Strand_t::StrandUnknown ||
         hit_strand == Strand_t::StrandUnknown;

         assert(iter_open->second.left_pos() == hit_partner_pos && strand_agree);
         if(iter_open->second._left_read == nullptr && iter_open->second._right_read){
            iter_open->second.set_left_read(move(hit));
         }
         else if(iter_open->second._right_read == nullptr && iter_open->second._left_read){
            iter_open->second.set_right_read(move(hit));
         }
         else{
            assert(false);
         }
         addHit(iter_open->second);
         _open_mates.erase(iter_open);
      }
   }
   return true;
}

void HitCluster::clearOpenMates()
{
   _open_mates.clear();
}

int HitCluster::collapseHits(){
   assert(_uniq_hits.empty());
   if(_hits.empty())
      return false;
   sort(_hits.begin(), _hits.end());
   _uniq_hits.resize(_hits.size());
   auto it_hits = _hits.begin();
   auto it_uniq = _uniq_hits.begin();
   while(it_hits != _hits.end()){
      *it_uniq = *it_hits;
      it_uniq->add_2_collapse_mass( it_uniq->raw_mass() );
      ++it_uniq;
      ++it_hits;
   }
   vector<PairedHit> duplicated;
   duplicated.reserve(_hits.size());
   auto last = unique2(_uniq_hits.begin(),_uniq_hits.end(), back_inserter(duplicated));
   _uniq_hits.erase(last, _uniq_hits.end());
   sort(duplicated.begin(), duplicated.end());
   size_t i = 0, j =0;
   while(i < _uniq_hits.size() && j < duplicated.size()){
      if(_uniq_hits[i] == duplicated[j]){
         _uniq_hits[i].add_2_collapse_mass( duplicated[j++].raw_mass() );
      }
      else
         ++i;
   }

//#ifdef DEBUG
//   if(duplicated.size() > 1){
//
//      for(auto &i: duplicated){
//         for(auto &j: _uniq_hits){
//            if(i == j) cout<<"left_pos "<<j.left_pos()<<"\t"<<j._collapse_mass<<endl;
//         }
//      }
//   }
//#endif
   assert(j == duplicated.size());
   return _uniq_hits.size();
}


bool HitCluster::overlaps( const HitCluster& rhs) const{
   if(ref_id() != rhs.ref_id()) return false;
   return left() < rhs.left() ? right() >= rhs.left() : left() <= rhs.right();
}

void HitCluster::guess_strand(){
   int plus_intron_size =0;
   int minus_intron_size =0;

     for(auto r = _uniq_hits.cbegin(); r< _uniq_hits.cend(); ++r){
        if(r->contains_splice()){
           if(r->strand() == Strand_t::StrandPlus)
              ++plus_intron_size;
           else
              ++minus_intron_size;
        }
     }
     if(plus_intron_size == 0 && minus_intron_size >0){
        _strand = Strand_t::StrandMinus;
     }
     else if(plus_intron_size > 0 && minus_intron_size ==0){
        _strand = Strand_t::StrandPlus;
     }
     else if(plus_intron_size == 0 && minus_intron_size == 0){
        _strand = Strand_t::StrandUnknown;
     }
     else{
        if(plus_intron_size / minus_intron_size > _kMinFold4BothStrand)
           _strand = Strand_t::StrandPlus;
        if(minus_intron_size / plus_intron_size > _kMinFold4BothStrand)
           _strand = Strand_t::StrandMinus;
        _strand = Strand_t::StrandBoth;
     }
}

bool ClusterFactory::loadRefmRNAs(vector<unique_ptr<GffSeqData>> &gseqs, RefSeqTable &rt,
      const char *seqFile)
{
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
            [](const unique_ptr<GffSeqData> &lhs, const unique_ptr<GffSeqData> &rhs){
               return lhs->get_gseq_id() < rhs->get_gseq_id();
      });
   }
   FaInterface fa_api;
   FaSeqGetter *fsg = NULL;

   // load sequencing if required
   if(seqFile != NULL){
      fa_api.initiate(seqFile);
   }
   for(uint i = 0; i<gseqs.size(); ++i){// for loop for each chromosome
      GffSeqData * gseq = &(*gseqs[i]);
      int f = 0;
      int r = 0;
      int u = 0;
      RefID ref_id = rt.get_id(gseqs[i]->_g_seq_name);
      GffmRNA *mrna = NULL;
      vector<Contig> ref_mrna_for_chr;
      if(fa_api.hasLoad()){
         delete fsg;
         fsg = NULL;
         fsg = new FaSeqGetter();
         fa_api.load2FaSeqGetter(*fsg,gseqs[i]->_g_seq_name);
         if(fsg == NULL){
            LOG_ERR("Reference sequence ", gseqs[i]->_g_seq_name, " can not be load!\n");
         }
      }
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
         Contig ref_contig(ref_id, strand, feats, true);
         ref_contig.annotated_trans_id(mrna->_transcript_id);
         //cout<<"ref contig left pos "<<ref_contig.left()<<endl;
         ref_mrna_for_chr.push_back(ref_contig);
      }// end while loop
      sort(ref_mrna_for_chr.begin(), ref_mrna_for_chr.end());
      _ref_mRNAs.insert(_ref_mRNAs.end(), ref_mrna_for_chr.begin(), ref_mrna_for_chr.end());
      ref_mrna_for_chr.clear();
   }//end for loop
   delete fsg;
   fsg = NULL;
   return true;
}

double ClusterFactory::next_valid_alignment(ReadHit& readin){
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
            const string cur_chr_name = _hit_factory->_ref_table.ref_name(readin.ref_id());
            const string last_chr_name = _hit_factory->_ref_table.ref_name(_prev_hit_ref_id);
            LOG_ERR("BAM file not sort correctly!");
            LOG_ERR("The current position is: ", cur_chr_name, ":", readin.left());
            LOG_ERR("and previous position is: ", last_chr_name, ":", _prev_hit_pos);
         }
      }

      _prev_hit_ref_id = readin.ref_id();
      _prev_hit_pos = readin.left();
      break;
   }
   return raw_mass;
}

double ClusterFactory::rewindHit(const ReadHit& rh)
{
   double mass = rh.mass();
   _hit_factory->undo_hit();
   return mass;
}

int ClusterFactory::addRef2Cluster(HitCluster &clusterOut){
   if(_refmRNA_offset >=  _ref_mRNAs.size()) {
      _has_load_all_refs = true;
      return 0;
   }

   clusterOut.addRefContig(&_ref_mRNAs[_refmRNA_offset++]);
   if(_refmRNA_offset >= _ref_mRNAs.size()){
      _has_load_all_refs = true;
      return 1;
   }

   size_t i = 0;
   while(i < clusterOut._ref_mRNAs.size()){
      Contig* ref = clusterOut._ref_mRNAs[i];
      if(Contig::overlaps_directional(*ref, _ref_mRNAs[_refmRNA_offset])){
         clusterOut.addRefContig(&_ref_mRNAs[_refmRNA_offset++]);
         if(_refmRNA_offset >= _ref_mRNAs.size()){
            _has_load_all_refs = true;
            return clusterOut._ref_mRNAs.size();
         }
         i=0;
      }
      else{
         ++i;
      }
   }
   return clusterOut._ref_mRNAs.size();
}

void ClusterFactory::rewindReference(HitCluster &clusterOut , int num_regress) {
   clusterOut.left(UINT_MAX);
   clusterOut.right(0);
   clusterOut.ref_id(-1);
   clusterOut._ref_mRNAs.clear();
   _refmRNA_offset -= num_regress;
   assert(_refmRNA_offset >= 0);
}

int ClusterFactory::nextCluster_denovo(HitCluster &clusterOut,
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

      if(clusterOut.size() == 0){ // add first hit
         clusterOut.addOpenHit(new_hit, true, true);
         clusterOut.addRawMass(mass);
//#ifdef DEBUG
//         cout<<clusterOut.right()<<"\t denovo cluster right"<<endl;
//#endif
      } else { //add the rest
         if(hit_lt_cluster(*new_hit, clusterOut, _kMaxOlapDist)){
            // should never reach here
            LOG_ERR("In alignments.cpp: It appears that SAM/BAM not sorted!");
         }
         if(hit_gt_cluster(*new_hit, clusterOut, _kMaxOlapDist)){
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


int ClusterFactory::nextCluster_refGuide(HitCluster &clusterOut)
{
   bool skip_read = false;
   if(!_hit_factory->recordsRemain()) return -1;
   // add first hit
   if(hasLoadRefmRNAs()){
      //if all ref mRNAs have been loaded but alignments haven't
      int num_added_refmRNA = addRef2Cluster(clusterOut);
      if( num_added_refmRNA == 0){
         return nextCluster_denovo(clusterOut);
      }
      //else
      // add as many as possible until encounter gap > _kMaxOlapDist
      while(true){
         ReadHitPtr new_hit(new ReadHit());
         double mass = next_valid_alignment(*new_hit);
         if(!_hit_factory->recordsRemain())
            return clusterOut.size();

         if(hit_lt_cluster(*new_hit, clusterOut, _kMaxOlapDist)){ // hit hasn't reach reference region
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
   //            cout<<"next:"<<next_ref_start_pos<<"\t hit chr:"<<new_hit->ref_id()<<"\thit left "<<new_hit->left()<<"\t previous:\t"<<clusterOut.left()<<endl;
   //#endif
               rewindReference(clusterOut, num_added_refmRNA);
               return nextCluster_denovo(clusterOut, next_ref_start_pos, next_ref_start_ref);
            }
         }

         if(hit_gt_cluster(*new_hit, clusterOut, _kMaxOlapDist)){ // read has gone to far.
            rewindHit(*new_hit);
            break;
         }

         clusterOut.addOpenHit(new_hit, false, false);
         clusterOut.addRawMass(mass);
         if(clusterOut._hits.size() >= HitCluster::_kMaxFragPerCluster ){
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dis(1, clusterOut._hits.size());
         }

      } // end while loop
   } // end loadRefmRNAs
   return clusterOut.size();
}

void ClusterFactory::mergeClusters(HitCluster & last, HitCluster &cur){
   /*
    * reassign paired hits;
    */
   sort(last._hits.begin(), last._hits.end());
   vector<PairedHit> imcompatible_hits;
   vector<vector<PairedHit>::iterator> imcomp_hits_its;
   bool first_incompatible = false;
   auto last_it = last._hits.begin();
   for(; last_it != last._hits.end(); ++last_it){
      if(last_it->contains_splice()){
         if(last_it->strand() != last.ref_strand()){
            first_incompatible = true;
            imcompatible_hits.push_back(*last_it);
            imcomp_hits_its.push_back(last_it);
         }
         else{
            first_incompatible = false;
         }
      }
      else{
         if(first_incompatible){
            imcompatible_hits.push_back(*last_it);
         }
      }
   }
   for(auto it = imcompatible_hits.cbegin(); it != imcompatible_hits.cend(); ++it){
      if(hit_complete_within_cluster(*it, cur, 0))
         cur._hits.push_back(*it);
   }
   for(auto &i: imcomp_hits_its){
      last._hits.erase(i);
   }

   /*
    * reassign open hits;
    */
   for(auto it = last._open_mates.begin(); it != last._open_mates.end(); ++it){
      if(it->second._left_read){
         cur.addOpenHit(it->second._left_read,false,false);
      }
      else{
         assert(it->second._right_read);
         cur.addOpenHit(it->second._right_read,false,false);
      }

   }
}

int ClusterFactory::ParseClusters(){
   int num_clusters = 0;
   vector<double> frag_len_dist;
   unique_ptr<HitCluster> last_cluster (new HitCluster());
   if( -1 == nextCluster_refGuide(*last_cluster) ) {
      return true;
   }
   while(true){
      unique_ptr<HitCluster> cur_cluster (new HitCluster());
      if(-1 == nextCluster_refGuide(*cur_cluster)){
         break;
      }
     if(cur_cluster->overlaps(*last_cluster) ){
        if(!cur_cluster->hasRefmRNAs() && !last_cluster->hasRefmRNAs()) {
           LOG_ERR("Error: It is unlikely that novo cluster overlaps");
           LOG_ERR(last_cluster->ref_id(), ":", last_cluster->left(), "-", last_cluster->right());
           LOG_ERR(cur_cluster->ref_id(),":", cur_cluster->left(), "-", cur_cluster->right());
        }
        mergeClusters(*last_cluster, *cur_cluster);
      }
     last_cluster->clearOpenMates();
     last_cluster->collapseHits();
     //last_cluster->guess_strand();
     size_t s = last_cluster->right() - last_cluster->left();
     vector<float> exon_doc(s, 0);
     if(last_cluster->size() > 0){
        compute_doc_4_cluster(*last_cluster, exon_doc);
     }
#ifdef DEBUG
     if(last_cluster->hasRefmRNAs()){
//           sort(last_cluster->_hits.begin(),last_cluster->_hits.end());
//           for(auto &h: last_cluster->_hits){
//              cout<<"spliced read on strand"<<h.strand()<<"\t"<< h.left_pos()<<"-"<<h.right_pos()<<endl;
//           }
            cout<<"number of Ref mRNAs "<<last_cluster->_ref_mRNAs.size()<<"\tRef cluster: "\
                  <<last_cluster->ref_id()<<"\t"<<last_cluster->left()<<"-"<<last_cluster->right()<<"\t"<<last_cluster->size()<<endl;
            cout<<"number of unique hits\t"<<last_cluster->_uniq_hits.size()<<endl;
         }
#endif
     last_cluster = move(cur_cluster);
   }
   last_cluster->clearOpenMates();
   last_cluster->collapseHits();
   //last_cluster->guess_strand();
}



void ClusterFactory::compute_doc_4_cluster(const HitCluster & hit_cluster, vector<float> &doc){

//#ifdef DEBUG
//   bool print  = false;
//   for(auto r = hit_cluster._uniq_hits.cbegin(); r< hit_cluster._uniq_hits.cend(); ++r){
//      for(auto &i: r->left_read_obj().cigar()){
//         if (i._type == INS)
//            print = true;
//      }
//      if(print){
//      Contig hit(*r);
//         cout<<"the read is "<<hit.left()<<"-"<<hit.right()<<endl;
//         cout<<"features are ";
//      for(auto i: hit._genomic_feats){
//         cout<<i._genomic_offset<<":"<<i._match_op._code<<":"<<i._match_op._len<<"; ";
//      }
//      cout<<endl;
//      print = false;
//      }
//   }
//#endif
     vector<Contig> hits;
     vector<IntronDOC> intron_counter;
     for(auto r = hit_cluster._uniq_hits.cbegin(); r< hit_cluster._uniq_hits.cend(); ++r){
        Contig hit(*r);
        hits.push_back(hit);
     }
     sort(hits.begin(), hits.end());
     compute_doc(hit_cluster.left(), hit_cluster.right(), hits, doc, intron_counter);
//     for(auto &i: doc){
//        cout<<"left: "<<i.left<<" right: "<<i.right<<" small: "<<i.small_span_read<<" total: "<<
//              i.total_junc_reads<<endl;
//
//     }
}

void ClusterFactory::compute_doc(const uint left, const uint right, const vector<Contig> & hits,
      vector<float> &exon_doc,  vector<IntronDOC> &intron_doc){
   assert(right > left);
   for(size_t i = 0; i<hits.size(); ++i){
      const vector<GenomicFeature> & g_feats = hits[i]._genomic_feats;
      for(size_t j = 0; j<g_feats.size(); ++j){
         const GenomicFeature & gf = g_feats[j];
         if( gf._match_op._code == Match_t::S_MATCH){
            for(size_t p = gf.left(); p < gf.right(); ++p){
               exon_doc[p-left] += hits[i].mass();
            }
         }
         else if( gf._match_op._code == Match_t::S_INTRON){
            IntronDOC cur_intron(gf.left(), gf.right());
            if(intron_doc.empty()){
               ++(cur_intron.total_junc_reads);
               if(g_feats[j-1]._match_op._len < kSmallOverHang &&
                     g_feats[j+1]._match_op._len < kSmallOverHang){
                  ++(cur_intron.small_span_read);
               }
               intron_doc.push_back(move(cur_intron));
               continue;
            }
            auto it = lower_bound(intron_doc.begin(), intron_doc.end(), cur_intron);
            if( *it == cur_intron){
               ++(it->small_span_read);
               if(g_feats[j-1]._match_op._len < kSmallOverHang &&
                     g_feats[j+1]._match_op._len < kSmallOverHang){
                  ++(it->small_span_read);
               }
            }
            else{
               ++(cur_intron.total_junc_reads);
               if(g_feats[j-1]._match_op._len < kSmallOverHang &&
                     g_feats[j+1]._match_op._len < kSmallOverHang){
                  ++(cur_intron.small_span_read);
               }
               intron_doc.insert(it, cur_intron);
            }
         }
      }
   }
}
