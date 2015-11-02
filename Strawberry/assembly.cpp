/*
 * assembly.cpp
 *
 *  Created on: Mar 18, 2015
 *      Author: ruolin
 */

#include "assembly.h"
#include "common.h"
#include "contig.h"
#include <iterator>
#include <iostream>
//#include <lemon/smart_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <lemon/core.h>

using namespace lemon;
typedef int LimitValueType;


bool FlowNetwork::comp_lt_first(const std::pair<uint, uint> & lhs, const std::pair<uint,uint> &rhs){
      return lhs.first < rhs.first;
   }

bool FlowNetwork::comp_lt_second(const std::pair<uint, uint> & lhs, const std::pair<uint,uint> &rhs){
      return lhs.second < rhs.second;
   }

bool FlowNetwork::search_left(const GenomicFeature &lhs, const uint rhs){
   return lhs._genomic_offset < rhs;
}

bool FlowNetwork::search_right(const GenomicFeature & lhs, const uint rhs){
   uint right = lhs._genomic_offset + lhs._match_op._len -1;
   return right < rhs;
}

//void FlowNetwork::initGraph(const int &left,
//        const vector<float> &exon_doc,
//        const vector<IntronTable> &intron_counter,
//        const vector<size_t> &bad_introns,
//        vector<GenomicFeature> &exons,
//        Graph::NodeMap<const GenomicFeature*> &node2feat)
//{
//   splicingGraph(left, exon_doc, intron_counter, bad_introns, exons);
//   createNetwork(exons, intron_counter, bad_introns, node2feat);
//}

// comparator used for searching in the split_bars object, which records the
// intron boundaries. We use pair<uint,bool> to indicate its position and
// whether it is left boundary or right boundary.
// left = true; right = false

void FlowNetwork::add_sink_source(Graph &g, Graph::Node &source, Graph::Node &sink){
   source = g.addNode();
   sink = g.addNode();
   for(Graph::NodeIt n(g); n != lemon::INVALID; ++n){
      if( n == source || n == sink )
         continue;
      int out_edges = 0;
      int in_edges = 0;
      for(Graph::OutArcIt out(g, n); out != lemon::INVALID; ++out)
         ++out_edges;
      for(Graph::InArcIt in(g,n); in != lemon::INVALID; ++in)
         ++in_edges;
      if(in_edges == 0){
         g.addArc(source, n);
      }
      if(out_edges == 0){
         g.addArc(n, sink);
      }
   }
   g.addArc(sink,source);
}



void FlowNetwork::flowDecompose(const Graph &g,
   const Graph::ArcMap<int> &flow,
   const Graph::Node &source,
   const Graph::Node &sink,
   std::vector<std::vector<Graph::Arc>> &paths ){

   Graph::ArcMap<int> copy_flow(g);
   for(Graph::ArcIt arc(g); arc != lemon::INVALID; ++arc){
      copy_flow[arc] = flow[arc];
   }

   while(hasFlow(g, copy_flow, source)){
      int bottle_neck = INT_MAX;
      std::vector<Graph::Arc> path;
      Graph::Node cur_node = source;
      while(cur_node != sink ){
         for(Graph::OutArcIt out(g, cur_node); out != lemon::INVALID; ++out){
            if(copy_flow[out] > 0){
               bottle_neck = max(bottle_neck, copy_flow[out]);
               cur_node = g.target(out);
               path.push_back(out);
               break;
            }
         }
      }
      for(auto edge: path){
         --copy_flow[edge];
      }
      paths.push_back(path);
   }
}


void FlowNetwork::splicingGraph(const int &left, const std::vector<float> &exon_doc,
      const std::map<pair<uint,uint>, IntronTable> &intron_counter,
      std::vector<GenomicFeature> &exons)
{
   /*
    * create non overlapping exon segments which will be
    * used as nodes in createNetwork();
    */

   vector<pair<uint,uint>> paired_bars; // two end intron boundaries
   vector<pair<uint,bool>> single_bars; // single intron boundaries

   for(auto i= intron_counter.cbegin(); i != intron_counter.cend(); ++i){
      paired_bars.push_back(pair<uint, uint> (i->first.first, i->first.second));
      single_bars.push_back(pair<uint, bool> (i->first.first, true));
      single_bars.push_back(pair<uint, bool> (i->first.second, false));
   }

   // unique element in single_bars
   sort(single_bars.begin(), single_bars.end(),comp_lt_first);
   auto newend = unique(single_bars.begin(), single_bars.end(),
         [](const pair<uint, bool> &lhs, const pair<uint, bool> &rhs){
         return lhs.first == rhs.first;}
         );
   single_bars.erase(newend, single_bars.end());

   // unique element in paired_bars
   sort(paired_bars.begin(), paired_bars.end());
   auto new_end = unique(paired_bars.begin(), paired_bars.end());
   paired_bars.erase(new_end, paired_bars.end());
   list<pair<uint,uint>> exon_boundaries;


   /*
    * preliminary exon segments.
    * */
   uint l = 0;
   uint r = l;
   for(size_t i = 0; i< exon_doc.size(); ++i){
      if(exon_doc[i] > 0 && l == 0){
         l = i+left;
      }
      if(exon_doc[i] == 0 && l != 0 ){
         r = i+left-1;
         exon_boundaries.push_back(pair<uint,uint>(l,r));
         l = 0;
      }
   }
   if( l != 0 && l < left+exon_doc.size() )
      exon_boundaries.push_back(pair<uint,uint>(l, left+exon_doc.size()-1));

   /*
    * When some exonic coverage gaps exist
    * due to low sequncing coverages.
    */
   auto iTer = exon_boundaries.begin();
   while(true){
      uint head = iTer->second;
      ++iTer;
      if(iTer == exon_boundaries.end()) break;
      uint tail = iTer->first;
      bool is_coverage_deficit = true;
      for(auto i= intron_counter.cbegin(); i != intron_counter.cend(); ++i){
         if(i->first.first-1 <= head && i->first.second+1 >= tail){
            is_coverage_deficit = false;
            break;
         }
      }
      if(is_coverage_deficit){
         iTer--;
         uint newStart = iTer->first;
         exon_boundaries.erase(iTer++);
         iTer->first = newStart;
      }
   };


   /*
    * for single exon genes
    * */

   if(paired_bars.size() == 0){
      uint l = exon_boundaries.front().first;
      uint r = exon_boundaries.back().second;
      exons.push_back(GenomicFeature(Match_t::S_MATCH, l, r-l+1));
      return;
   }

   /*
    * further divided preliminary exon segments into smaller pieces based on intorn boundaries.
    * */
   auto e = exon_boundaries.begin();
   size_t s = 0;
   while(e != exon_boundaries.end() && s < single_bars.size()){

      uint bar = single_bars[s].first;
      bool left = single_bars[s].second;
      if(bar < e->first){
         ++s;
      }
      else if(bar >= e->first && bar <= e->second){
         uint temp = e->second;
         if(left){
            e->second = bar-1;
            exon_boundaries.insert(++e, pair<uint,uint>(bar, temp));
         }
         else{
            e->second = bar;
            exon_boundaries.insert(++e, pair<uint,uint>(bar+1, temp));
         }
         --e;
         ++s;
      }
      else{
         ++e;
      }
   }

   /*
   * filter exon segments if it does not have intron supporting
   * */
   vector<size_t> dropoff;
   vector<pair<uint, uint>> left_coords;
   vector<pair<uint, uint>> right_coords;
   vector<pair<uint, uint>> e_boundaries;
   for(auto i: exon_boundaries){
      e_boundaries.push_back(i);
   }

   for(uint i = 0; i < paired_bars.size(); ++i){
      left_coords.push_back(pair<uint, uint>(paired_bars[i].first, i));
      right_coords.push_back(pair<uint, uint>(paired_bars[i].second, i));
   }
   sort(left_coords.begin(), left_coords.end(), comp_lt_first);
   sort(right_coords.begin(), right_coords.end(), comp_lt_first);

   for(size_t ex = 0 ; ex < e_boundaries.size(); ++ex){
      if(ex == 0){
         if(e_boundaries[ex].second == e_boundaries[ex+1].first-1);
         else{
            auto lower = lower_bound(left_coords.begin(), left_coords.end(), pair<uint, uint>(e_boundaries[ex].second+1,0), comp_lt_first);
            if(lower != left_coords.end() && lower->first == e_boundaries[ex].second+1){
               uint right = paired_bars[lower->second].second;
               if(!binary_search(e_boundaries.begin(), e_boundaries.end(), pair<uint, uint>(right+1,0), comp_lt_first)){
                  dropoff.push_back(ex);
               }
            }
            else{
               dropoff.push_back(ex);
            }
         }
      }
      else if(ex == e_boundaries.size()-1){
         if( e_boundaries[ex].first == e_boundaries[ex-1].second+1 );
         else{
            auto lower = lower_bound(right_coords.begin(), right_coords.end(), pair<uint,uint>(e_boundaries[ex].first-1,0), comp_lt_first);
            if( lower != right_coords.end() && lower->first == e_boundaries[ex].first-1){
               uint left = paired_bars[lower->second].first;
               if(!binary_search(e_boundaries.begin(), e_boundaries.end(), pair<uint, uint>(0,left-1), comp_lt_second))
                  dropoff.push_back(ex);
            }
            else
               dropoff.push_back(ex);
         }
      }
      else{
         if( e_boundaries[ex].first != e_boundaries[ex-1].second+1 || e_boundaries[ex].second != e_boundaries[ex+1].first-1){
            bool no_intron_on_left = false;
            bool no_intron_on_right = false;
            auto l = lower_bound(left_coords.begin(), left_coords.end(), pair<uint, uint>(e_boundaries[ex].second+1,0), comp_lt_first);
            if(l != left_coords.end() && l->first == e_boundaries[ex].second+1){
               uint right = paired_bars[l->second].second;
               if(!binary_search(e_boundaries.begin(), e_boundaries.end(), pair<uint, uint>(right+1,0), comp_lt_first))
                  no_intron_on_right = true;

            }
            else{
               no_intron_on_right = true;
            }

            auto r = lower_bound(right_coords.begin(), right_coords.end(), pair<uint,uint>(e_boundaries[ex].first-1,0), comp_lt_first);
            if( r != right_coords.end() && r->first == e_boundaries[ex].first-1){
               uint left = paired_bars[r->second].first;
               if(!binary_search(e_boundaries.begin(), e_boundaries.end(), pair<uint, uint>(0,left-1), comp_lt_second))
                  no_intron_on_left = true;
            }
            else{
               no_intron_on_left = true;
            }
            if(no_intron_on_left && no_intron_on_right){
               dropoff.push_back(ex);
            }
         }
      } // end of if-ifelse-else condition
//      else{
//
//         // first two conditiosn cover situations when exon
//         // segments are supported by intron only on one side.
//         //
//         if( next(e)->first - e->second > 1 &&
//             !binary_search(paired_bars.begin(), paired_bars.end(), pair<uint, bool>(e->second+1,true), comp_lt)
//           )
//         {
//            if( binary_search(paired_bars.begin(), paired_bars.end(), pair<uint, bool>(e->first-1,true), comp_lt)){
//               if(e->first - prev(e)->second == 2) dropoff.push_back(e);
//            }
//            else{
//               dropoff.push_back(e);
//            }
//            continue;
//         }
//
//         if (e->first - prev(e)->second > 1 &&
//               !binary_search(paired_bars.begin(), paired_bars.end(), pair<uint, bool>(e->first-1,true), comp_lt)
//            )
//         {
//            if(binary_search(paired_bars.begin(), paired_bars.end(), pair<uint, bool>(e->second+1,true), comp_lt)){
//               if(next(e)->first - e->second == 2) dropoff.push_back(e);
//            }
//            else{
//               dropoff.push_back(e);
//            }
//            continue;
//         }
//
//         if( ( binary_search(paired_bars.begin(), paired_bars.end(), pair<uint, bool>(e->first-1,true), comp_lt) ||
//               e->first == prev(e)->second+1
//             ) &&
//             (
//               binary_search(paired_bars.begin(), paired_bars.end(), pair<uint, bool>(e->second+1,true), comp_lt) ||
//               e->second == next(e)->first-1
//             ));
//         else{
//            dropoff.push_back(e);
//         }
//      }

   }

//#ifdef DEBUG
//   for(auto i: exon_boundaries)
//      cout<<"left: "<<i.first<<" right: "<<i.second<<endl;
//   for(auto s: paired_bars)
//      cout<<"intorn: "<<s.first<<"-"<<s.second<<endl;
//   cout<<"---------------------"<<endl;
//#endif

   vector<list<pair<uint,uint>>::iterator> drops;
   for(auto d: dropoff){
      drops.push_back(next(exon_boundaries.begin(),d));
   }

   for(auto&d: drops){
      exon_boundaries.erase(d);
   }


   for(auto i: exon_boundaries){
      exons.push_back(GenomicFeature(Match_t::S_MATCH, i.first, i.second-i.first+1));
      double sum = accumulate(exon_doc.begin() + i.first - left, exon_doc.begin() + i.second+1-left,0);
      double avg_doc = sum/(i.second - i.first +1);
      exons.back().avg_doc(avg_doc);
   }
   sort(exons.begin(), exons.end());
}

bool FlowNetwork::createNetwork(
      const vector<Contig> &hits,
      const vector<GenomicFeature> &exons,
      const std::map<pair<uint,uint>, IntronTable> &intron_counter,
      const vector<vector<size_t>> &constraints,
      Graph::NodeMap<const GenomicFeature*> &node2feat,
      Graph::ArcMap<int> &cost_map,
      Graph::ArcMap<int> &min_flow_map,
      vector<vector<Graph::Arc>> &path_cstrs)
{
   assert(!hits.empty());
   vector<Graph::Node> nodes;
   vector<Graph::Arc> arcs;
   map<const GenomicFeature*, Graph::Node> feat2node;
   for(size_t i = 0; i< exons.size(); ++i){
      nodes.push_back(_g.addNode());
      node2feat[nodes[i]] = &exons[i];
      feat2node[&exons[i]] = nodes[i];
   }

   // single exon case;
   if(exons.size() == 1) return true;

   // 1) add arc defined by intron
   for(auto i = intron_counter.cbegin(); i != intron_counter.cend(); ++i){
      const IntronTable &intron = i->second;
      auto e1 = lower_bound(exons.begin(), exons.end(), intron.left-1, search_right);
      auto e2 = lower_bound(exons.begin(), exons.end(), intron.right+1, search_left);
      if(e1 == exons.end() || e2 == exons.end()) {
#ifdef DEBUG
         assert(false);
#endif
         return false;
      }

      //cout<< intron.left-1 << " vs " << intron.right+1<<endl;
      //cout<< e1->right()<<" and "<< e2->left()<<endl;
      Graph::Arc a = _g.addArc( feat2node[&(*e1)],  feat2node[&(*e2)]);
      arcs.push_back(a);
   }

   // 2) add arc defined by consecutive exons
   for(size_t i=0; i < exons.size()-1; ++i){
      uint right = exons[i]._genomic_offset + exons[i]._match_op._len -1;
      if(exons[i+1]._genomic_offset - right == 1){
         const Graph::Node &node_a = feat2node[&exons[i]];
         const Graph::Node &node_b = feat2node[&exons[i+1]];
         Graph::Arc a = _g.addArc(node_a, node_b);
         arcs.push_back(a);
      }
   }
   addWeight(hits, intron_counter, node2feat, cost_map);

   //Now add subpath constraints to the just created graph.
   InDegMap<Graph> inDeg(_g);
   OutDegMap<Graph> outDeg(_g);

   for(auto c: constraints){
      vector<Graph::Arc> path_cstr;

      size_t s_idx = c.front();
      size_t t_idx = c.back();
      const Graph::Node &s = feat2node[&exons[s_idx]];
      const Graph::Node &t = feat2node[&exons[t_idx]];

      //At least one of the internal nodes have >1 out-degree and >1 in-degree
      bool is_valid = false;
      for(size_t idx = 1; idx<c.size()-1; ++idx){
         const Graph::Node &n = feat2node[&exons[c[idx]]];
//         if(exons[c[idx]]._genomic_offset == 7384){
//            cout<<"haha"<<inDeg[n]<<"\t"<<outDeg[n]<<endl;
//         }

         if(inDeg[n] > 1 && outDeg[n] > 1)
            is_valid = true;
      }

      if(ArcLookUp<ListDigraph>(_g)(s,t) == INVALID && is_valid){

         LOG("Path Constraints:\t");
         for(size_t i = 0; i< c.size()-1; ++i){
            const Graph::Node &pre = feat2node[&exons[c[i]]];
            const Graph::Node &sec = feat2node[&exons[c[i+1]]];
            LOG("(", exons[c[i]].left(), ",",exons[c[i]].right(),\
                ")-(", exons[c[i+1]].left(),",",exons[c[i+1]].right(),")\t");
//#ifdef DEBUG
//            cout<<exons[c[i]]._genomic_offset<<"---"<<exons[c[i+1]]._genomic_offset<<endl;
//#endif
            Graph::Arc arc_found = ArcLookUp<ListDigraph>(_g)(pre, sec);
            if(arc_found == INVALID){

               LOG_ERR("Calculating Path Constraints failed");
               LOG_ERR("No intron connects exon ", hits[0].ref_id(), ":",node2feat[pre]->left(), "-", node2feat[pre]->right());
               LOG_ERR("No intron connects exon ", hits[0].ref_id(), ":",node2feat[sec]->left(), "-", node2feat[sec]->right());
               continue;
            }

            path_cstr.push_back(arc_found);
         }
         path_cstrs.push_back(path_cstr);
      }
   }

   if(path_cstrs.empty()){
      for(auto arc: arcs){
         min_flow_map[arc] = 1;
      }
      return true;
   }

   vector<Graph::Arc> one_d_path_cstrs;
   for(auto p: path_cstrs){
     for(auto e: p)
        one_d_path_cstrs.push_back(e);
   }
   sort(one_d_path_cstrs.begin(), one_d_path_cstrs.end());
   auto new_end = unique(one_d_path_cstrs.begin(), one_d_path_cstrs.end());
   one_d_path_cstrs.erase(new_end, one_d_path_cstrs.end());

   for(auto edge: arcs){
      if(find(one_d_path_cstrs.begin(), one_d_path_cstrs.end(), edge) == one_d_path_cstrs.end()){
         vector<Graph::Arc> path;
         path.push_back(edge);
         path_cstrs.push_back(path);
      }
   }

   for(auto p: path_cstrs){
      assert(!p.empty());
      if(p.size()>1){
         // a subpath
         int cost = 0;
         for(auto arc : p){
            cost += cost_map[arc];
         }
         const Graph::Node s = _g.source(p.front());
         const Graph::Node t = _g.target(p.back());
         if( findArc(_g, s,t) == INVALID){
            const Graph::Arc a = _g.addArc(s,t);
            cost_map[a] = cost;
            min_flow_map[a] = 1;
         }
      }
      else{
         // or a arc
         min_flow_map[p[0]] = 1;
      }
   }
   return true;
}

void FlowNetwork::addWeight(const vector<Contig> &hits,
      const std::map<pair<uint,uint>, IntronTable> &intron_counter,
      const Graph::NodeMap<const GenomicFeature*> &node_map,
      Graph::ArcMap<int> &arc_map)
{
   // the following has O(n^2) complicity. We may need better solution in the future.
   for(Graph::ArcIt arc(_g); arc != lemon::INVALID; ++arc){
      const Graph::Node &s = _g.source(arc);
      const Graph::Node &t = _g.target(arc);
      uint arc_s = node_map[s]->right();
      uint arc_e = node_map[t]->left();
      float num_read_support = 0;
      if(arc_e - arc_s == 1){
         for(auto mp: hits){
            if(mp.left() > arc_e) break;
            if(mp.right() < arc_s) continue;
            for(auto feature: mp._genomic_feats){
               if(feature._match_op._code == Match_t::S_MATCH){
                  if(feature.left() <= arc_s-kMinDist4ExonEdge &&
                     feature.right() >= arc_e+kMinDist4ExonEdge){
                     num_read_support += mp.mass();
                  }
               }
            }
         }
      }
      else{
         arc_s += 1;
         arc_e -= 1;
         for(auto i  = intron_counter.cbegin(); i != intron_counter.cend(); ++i){
            if(arc_s == i->first.first && arc_e == i->first.second){
               num_read_support = i->second.total_junc_reads;
               break;
            }
         }
      }
      _max_weight = max(_max_weight, num_read_support);
      arc_map[arc] = num_read_support;
   }

   for(Graph::ArcIt arc(_g); arc != lemon::INVALID; ++arc){
      arc_map[arc] = _max_weight - arc_map[arc];
      assert(arc_map[arc] >= 0.0);
   }
}

void FlowNetwork::assignExonBin(
      const vector<GenomicFeature> & exons,
      const vector<Contig> &hits,
      const vector<Isoform> &transcripts,
      map<set<uint>, ExonBin> & exon_bin_map)
{
   for(auto mp = hits.cbegin(); mp != hits.cend(); ++mp){
      ExonBin eb(mp->ref_id());

      // first loop finds exon bin
      for(size_t i = 0; i< exons.size(); ++i){
         if(Contig::overlaps_only_on_exons(*mp, exons[i])){
            bool status = eb.insert_exon(&exons[i]);
            if(!status){
               LOG_ERR("Error in FlowNetwork::assignExonBin. A read overlaps a exon twice!!");
            }
         }
      } // end inner for loop
      if(eb._exon_in_bin.size() == 0){
         continue;
//#ifdef DEBUG
//         for(auto &e : exons){
//            cout<<"exon: "<<e.left()<<"-"<<e.right()<<endl;
//         }
//         cout<<"hit: chr: "<<mp->ref_id()<<"\t"<<mp->left()<<"-"<<mp->right()<<endl;
//#endif
      }
      eb.add_hit(&(*mp));

      map<set<uint>, ExonBin>::iterator it_ret = exon_bin_map.find(eb._coords);
      if(it_ret == exon_bin_map.end()){
         eb.add_isoform(transcripts);
         bool status;
         map<set<uint>, ExonBin>::iterator it;
         tie(it, status) = exon_bin_map.insert(pair<set<uint>, ExonBin> (eb._coords, move(eb)));
         assert(status);
      }
      else{
         it_ret->second.read_num_increase_by_1();
         it_ret->second.add_hit(eb);
      }
   }


   for(auto kv = exon_bin_map.cbegin(); kv != exon_bin_map.cend(); ++kv){

   }
}


vector<vector<size_t>> FlowNetwork::findConstraints(

   const vector<GenomicFeature> &exons,
   const vector<Contig> &hits)
{
   vector<vector<size_t>> result;
   for(auto mp = hits.cbegin(); mp != hits.cend(); ++mp){
      vector<size_t> constraint;
      vector<GenomicFeature> introns;
      uint start = 0;
      uint end = 0;
      for(auto feat = mp->_genomic_feats.cbegin(); feat != mp->_genomic_feats.cend(); ++feat){
         if(feat->_match_op._code != Match_t::S_GAP){
            if(feat->_match_op._code == Match_t::S_INTRON){
               introns.push_back(*feat);
            }
            if(start == 0){
               start = feat->left();
               end = feat->right();
            }
            else{
               end = feat->right();
            }
         }
         // if Match_t::S_GAP
         else{
            assert(end>start);
            for(size_t i = 0; i< exons.size(); ++i){
               if(GenomicFeature::overlap_in_genome(exons[i], start, end)){
                  bool valid = true;
                  for(auto intron: introns){
                     if(intron.contains(exons[i])){
                        valid = false;
                        break;
                     }
                  }

                  if(valid)
                     constraint.push_back(i);

               }
            }
            if(constraint.size() > 2){
               result.push_back(constraint);
            }
            start = 0;
            constraint.clear();
            introns.clear();
         }
         // repeat Match_t::S_GAP
         if(*feat == mp->_genomic_feats.back()){
            assert(end>start);
            for(size_t i = 0; i< exons.size(); ++i){
               if(GenomicFeature::overlap_in_genome(exons[i], start, end)){
                  bool valid = true;
                  for(auto intron: introns){
                     if(intron.contains(exons[i])){

                        valid = false;
                        break;
                     }
                  }
                  if(valid)
                     constraint.push_back(i);
               }
            }
            if(constraint.size() > 2){
               result.push_back(constraint);
            }
         }
      }
   }
   sort(result.begin(), result.end());
   auto new_end = unique(result.begin(), result.end());
   result.erase(new_end, result.end());
   return result;
}

bool FlowNetwork::solveNetwork(const Graph::NodeMap<const GenomicFeature*> &node_map,
      const vector<GenomicFeature> &exons,
      const vector<vector<Graph::Arc>> &path_cstrs,
      Graph::ArcMap<int> &cost_map,
      Graph::ArcMap<int> &min_flow_map,
      vector<vector<GenomicFeature>> &transcripts){

   // single exon cases
   if(exons.size() == 1){
      transcripts.push_back(exons);
   }

   add_sink_source(_g, _source, _sink);
   Graph::ArcMap<LimitValueType> u(_g);
   NetworkSimplex<Graph, LimitValueType, LimitValueType> FlowNetwork(_g);

   for(Graph::ArcIt arc(_g); arc != INVALID; ++arc){
      u[arc] = FlowNetwork.INF;
   }
   FlowNetwork.lowerMap(min_flow_map).upperMap(u).costMap(cost_map);

   NetworkSimplex<Graph>::ProblemType ret = FlowNetwork.run();

   // get flow of arcs
   Graph::ArcMap<LimitValueType> flow(_g);
   FlowNetwork.flowMap(flow);
   if(ret == NetworkSimplex<Graph>::INFEASIBLE || ret == NetworkSimplex<Graph>::UNBOUNDED){
#ifdef DEBUG
      digraphWriter(_g).                  // write g to the standard output
        arcMap("cost", cost_map).          // write 'cost' for for arcs
        arcMap("flow", min_flow_map).          // write 'flow' for for arcs
        node("source", _source).             // write s to 'source'
        node("target", _sink).             // write t to 'target'
        run();
      fprintf(stderr, "Infeasible or unbounded FlowNetwork flow\n");
      assert(false);
#endif
      return false;
   }
   vector<vector<Graph::Arc>> paths;
   flowDecompose(_g, flow, _source, _sink, paths);

   //put paths into vector<vector<GenomicFeature>>
   for(auto p: paths){
      vector<GenomicFeature> tscp;
      for(size_t i = 1; i <p.size(); ++i){
         const Graph::Arc e = p[i];
         const Graph::Node &arc_s = _g.source(e);
         const Graph::Node &arc_t = _g.target(e);
         bool is_edge = true;
         // this loop transform the path constrains to original meaning: a set of arcs.
         for(auto cstr: path_cstrs){
            const Graph::Node &path_s = _g.source(cstr.front());
            const Graph::Node &path_t = _g.target(cstr.back());
            if(arc_s == path_s && arc_t == path_t){

               is_edge = false;

               for(size_t idx = 0; idx<cstr.size()-1; ++idx){
                  const Graph::Node &n1 = _g.source(cstr[idx]);
                  const Graph::Node &n2 = _g.source(cstr[idx+1]);
                  tscp.push_back(GenomicFeature(Match_t::S_MATCH, node_map[n1]->_genomic_offset, node_map[n1]->_match_op._len));
                  if(idx+1 < cstr.size()-1)
                     if(node_map[n2]->left()-node_map[n1]->right() > 1){
                        tscp.push_back(GenomicFeature(Match_t::S_INTRON, node_map[n1]->right()+1, node_map[n2]->left()-1-node_map[n1]->right()));
                     }
               }

               const Graph::Node &n1 = _g.source(cstr.back());
               const Graph::Node &n2 = _g.target(cstr.back());
               tscp.push_back(GenomicFeature(Match_t::S_MATCH, node_map[n1]->_genomic_offset, node_map[n1]->_match_op._len));
               if(node_map[n2]->left()-node_map[n1]->right() > 1){
                  tscp.push_back(GenomicFeature(Match_t::S_INTRON, node_map[n1]->right()+1, node_map[n2]->left()-1-node_map[n1]->right()));
               }

               break;
            }
         }
         // else it is edge originally
         if(is_edge){
            //cout<<node_map[s]->left()<<"-"<<node_map[s]->right()<<"\t";
            tscp.push_back(GenomicFeature(Match_t::S_MATCH, node_map[arc_s]->_genomic_offset, node_map[arc_s]->_match_op._len));
            if( i+1 < p.size()){
               if(node_map[arc_t]->left()-node_map[arc_s]->right() > 1){
                  tscp.push_back(GenomicFeature(Match_t::S_INTRON, node_map[arc_s]->right()+1, node_map[arc_t]->left()-1-node_map[arc_s]->right()));
               }
            }
         }
      }
      transcripts.push_back(tscp);
   }
   return true;
}
