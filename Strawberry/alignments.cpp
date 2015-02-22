/*
 * alignment.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: ruolin
 */

#include "alignments.h"
#include "fasta.h"
#include <algorithm>
#include <random>
#ifdef DEBUG
   #include <iostream>
#endif
using namespace std;

int HitCluster::_next_id =0;
HitCluster::HitCluster(): _id(++_next_id){}

RefID HitCluster::ref_id() const{
   return _ref_id;
}

void HitCluster::addRefContig(Contig * contig){
   if(_ref_id != -1){
      assert(_ref_id == contig->ref_id());
   }
   else{
      _ref_id = contig->ref_id();
   }
   _leftmost = min(_leftmost, contig->left());
   _rightmost = max(_rightmost, contig->right());
   _ref_contigs.push_back(contig);
}

uint HitCluster::left() const {
   return _leftmost;
}
uint HitCluster::right() const{
   return _rightmost;
}

int HitCluster::size() const {
   return _hits.size();
}

bool HitCluster::addHit(const shared_ptr<PairedHit> &hit){
   if(_final){
      return false;
   }
   _leftmost = min(_leftmost, hit->left_pos());
   _rightmost = max(_rightmost, hit->right_pos());
   _hits.push_back(hit);
   return true;
}


bool HitCluster::addOpenHit(shared_ptr<ReadHit> hit)
{
   if(hit->right() - hit->left() > _kMaxGeneLen){
      SMessage("Warning: hit start at %d is longer than max gene length, skipping\n", hit->left());
      return false;
   }
   if(hit->is_singleton()){
      _ref_id = hit->interval().seq_id();
      shared_ptr<PairedHit> ph( new PairedHit(hit, nullptr));
      addHit(ph);
   }
   else{
     unordered_map<int, shared_ptr<PairedHit>>::iterator iter_open = _open_mates.find(hit->left());
      if( iter_open == _open_mates.end()){
         if(hit->left() <= hit->partner_pos()){
            shared_ptr<PairedHit> open_hit (new PairedHit(hit, nullptr));
            unordered_map<int, shared_ptr<PairedHit>>::iterator ins_pos;
            bool status;
            tie(ins_pos, status) = _open_mates.insert(make_pair(hit->partner_pos(), open_hit));
            if(!status) SMessage("Warning: Inserting read into open_mats failed\n");
         } else{
            return false;
         }
      } else{
         bool found_mate = false;
         if(iter_open->second->ref_seq_id() != hit->ref_id())
            return false;

         bool strand_agree = iter_open->second->strand() == hit->strand() ||
               iter_open->second->strand() == GenomicInterval::kStrandUnknown ||
               hit->strand() == GenomicInterval::kStrandUnknown;
         if(iter_open->second->left_pos() == hit->partner_pos() && strand_agree){
            iter_open->second->set_right_read(move(hit));
            addHit(iter_open->second);
            _open_mates.erase(iter_open);
            found_mate = true;
         }

         if(!found_mate){
            if(hit->left() <= hit->partner_pos()){
               shared_ptr<PairedHit> open_hit (new PairedHit(hit, nullptr));
              unordered_map<int, shared_ptr<PairedHit>>::iterator ins_pos;
               bool status;
               tie(ins_pos, status) = _open_mates.insert(make_pair(hit->partner_pos(), open_hit));
               if(!status) SMessage("Warning: Inserting read into open_mats failed\n");
            } else{
               return false;
            }
         }// end if(!found_partern)
      }

   }
   return true;
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
            SMessage("Warning: Sam file and Gff file are not sorted in the same order!\n");
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
            SMessage("Reference sequence %s can not be load!\n",gseqs[i]->_g_seq_name.c_str());
         }
      }
      int f_total = gseqs[i]->_forward_rnas.size();
      int r_total = gseqs[i]->_reverse_rnas.size();
      int u_total = gseqs[i]->_unstranded_rnas.size();
      while(!(f==f_total && r==r_total && u == u_total)){
         char strand = 0;
         if(f <f_total){
            mrna = &(*gseqs[i]->_forward_rnas[f++]);
            strand = '+';
         } else if(r < r_total){
            mrna = &(*gseqs[i]->_reverse_rnas[r++]);
            strand = '-';
         } else{
            mrna = &(*gseqs[i]->_unstranded_rnas[u++]);
            strand = '.';
         }

         vector<GenomicFeature> feats;
         for(uint e = 0; e < mrna->_exons.size(); ++e){
            GffExon& ex = *(mrna->_exons[e]);
            feats.push_back(GenomicFeature(S_MATCH, ex._iv.left(), ex._iv.right()-ex._iv.left()+1));
            if( e + 1 < mrna->_exons.size()){
               GffExon& next_ex  = *(mrna->_exons[e+1]);
               feats.push_back(GenomicFeature(S_INTRON, ex._iv.right()+1, next_ex._iv.left()-1-ex._iv.right() ));
            }
         }
         Contig ref_contig(ref_id, strand, feats, true);
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
            SError("Error:BAM file not sort correctly! The current position is %s:%d and previous position is %s:%d.\n",
                  cur_chr_name.c_str(), readin.left(), last_chr_name.c_str(), _prev_hit_pos);
         }
      }

      _prev_hit_ref_id = readin.ref_id();
      _prev_hit_pos = readin.left();
      break;
   }
   return raw_mass;
}

double ClusterFactory::rewind_hit(const ReadHit& rh)
{
   double mass = rh.mass();
   _hit_factory->undo_hit();
   return mass;
}

int ClusterFactory::addRef2Cluster(HitCluster &clusterOut){
   if(_refmRNA_offset >= _ref_mRNAs.size())
      return 0;
   clusterOut.addRefContig(&_ref_mRNAs[_refmRNA_offset++]);
   for(auto ref: clusterOut._ref_contigs){
      if(Contig::overlaps(*ref, _ref_mRNAs[_refmRNA_offset])){
         clusterOut.addRefContig(&_ref_mRNAs[_refmRNA_offset++]);
      }
   }
   return clusterOut._ref_contigs.size();
}

int ClusterFactory::nextCluster_denovo(HitCluster &clusterOut){
   if(!_hit_factory->recordsRemain()) return -1;
   while(true){
      shared_ptr<ReadHit> new_hit(new ReadHit());
      double mass = next_valid_alignment(*new_hit);
      //cout<<new_hit->left()<<endl;
      if(!_hit_factory->recordsRemain()){
         return clusterOut.size();
      }
      if(clusterOut.size() == 0){ // add first hit
         clusterOut.addOpenHit(new_hit);
      } else { //add the rest
         if(hit_lt_cluster(*new_hit, clusterOut)){
            // read hasn't reach reference region
            continue;
         }
         if(hit_gt_cluster(*new_hit, clusterOut)){
            // read has gone to far.
            rewind_hit(*new_hit);
            break;
         }
         clusterOut.addOpenHit(new_hit);
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
      // all ref mRNAs have been loaded but alignments haven't
      if( addRef2Cluster(clusterOut) == 0){
         //return -1;
         return nextCluster_denovo(clusterOut);
      }
   // add as many as possible until encounter gap > _kMinOlapDist
      while(true){
         shared_ptr<ReadHit> new_hit(new ReadHit());
         double mass = next_valid_alignment(*new_hit);
         if(!_hit_factory->recordsRemain())
            return clusterOut.size();
         if(hit_lt_cluster(*new_hit, clusterOut)){
            // read hasn't reach reference region
            continue;
         }
         if(hit_gt_cluster(*new_hit, clusterOut)){
            // read has gone to far.
            rewind_hit(*new_hit);
            break;
         }
         clusterOut.addOpenHit(new_hit);
         if(clusterOut._hits.size() >= HitCluster::_kMaxFragsSize ){
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dis(1, clusterOut._hits.size());
         }

      } // end while loop
   } // end loadRefmRNAs
   return clusterOut.size();
}


bool hit_lt_cluster(const ReadHit& hit, const HitCluster& cluster){
   if(hit.ref_id() != cluster.ref_id())
      return hit.ref_id() < cluster.ref_id();
   else
      return hit.right() < cluster.left();
}

bool hit_gt_cluster(const ReadHit& hit, const HitCluster& cluster){
   if(hit.ref_id() != cluster.ref_id()){
      //cout<<"shouldn't\t"<<hit.ref_id()<<":"<<cluster.ref_id()<<endl;
      return hit.ref_id() > cluster.ref_id();
   }
   else{
      return hit.left() > cluster.right();
   }
}
