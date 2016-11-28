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

#include "assembly.h"
#include "contig.h"
#include <iterator>
#include <iostream>
#include <climits>
//#include <lemon/smart_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <lemon/dijkstra.h>
#include <lemon/core.h>
#include <stack>

using namespace lemon;
typedef int LimitValueType;

void compute_exon_doc(const int left, const std::vector<float> & exon_doc, std::vector<GenomicFeature>& exons)
{
   for(uint i = 0; i != exons.size(); ++i){
      auto it_start = next(exon_doc.begin(), exons[i].left() - left);
      auto it_end = next(exon_doc.begin(), exons[i].right() - left);
      double cov = accumulate(it_start, it_end, 0.0);
      exons[i].avg_doc( cov/exons[i].len() );
   }
}


void  constraints_4_single_end(const std::vector<GenomicFeature>& feats,
                               const std::vector<GenomicFeature> & exons,
                               std::vector<size_t> &constraint
                              )
{
   uint start = 0;
   uint end = 0;
   std::vector<GenomicFeature> introns;
   for(auto feat = feats.cbegin(); feat != feats.cend(); ++feat){
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
}


bool FlowNetwork::comp_lt_first(const std::pair<uint, uint> & lhs, const std::pair<uint,uint> &rhs)
{
   return lhs.first < rhs.first;
}

bool FlowNetwork::comp_lt_second(const std::pair<uint, uint> & lhs, const std::pair<uint,uint> &rhs)
{
   return lhs.second < rhs.second;
}

bool FlowNetwork::search_left(const GenomicFeature &lhs, const uint rhs)
{
   return lhs._genomic_offset < rhs;
}

bool FlowNetwork::search_right(const GenomicFeature & lhs, const uint rhs)
{
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

void FlowNetwork::add_sink_source(Graph &g, Graph::Node &source, Graph::Node &sink)
{
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
   const Graph::ArcMap<int> &cost,
   const Graph::Node &source,
   const Graph::Node &sink,
   std::vector<std::vector<Graph::Arc>> &paths )
{

   Graph::ArcMap<int> copy_flow(g);
   Graph::ArcMap<int> edge_cost(g);
   for(Graph::ArcIt arc(g); arc != lemon::INVALID; ++arc){
      copy_flow[arc] = flow[arc];
   }
   for(Graph::ArcIt arc(g); arc != lemon::INVALID; ++arc){
      edge_cost[arc] = cost[arc];
   }
   /*assign weight to source node edges*/
   for(Graph::OutArcIt out(g, source); out != lemon::INVALID; ++out){
      int opt_cost = INT_MAX;
      Graph::Node cur_node = g.target(out);
      for(Graph::OutArcIt out2(g, cur_node); out2 != lemon::INVALID; ++out2){
         opt_cost = std::min(opt_cost, cost[out2]);
      }
      edge_cost[out] = opt_cost;
   }

   /*Greed way to decompose flow
    * In the future, maybe consider dynamic algorithm. But need to
    * average edge costs. If not, the path with fewer edges will be
    * chosen.
    * */
   while(hasFlow(g, copy_flow, source)){
//      int bottle_neck = INT_MAX;
      std::vector<Graph::Arc> path;
      Graph::Node cur_node = source;
      while(cur_node != sink ){
         Graph::Arc opt_arc = lemon::INVALID;
         int opt_cost = INT_MAX;
         for(Graph::OutArcIt out(g, cur_node); out != lemon::INVALID; ++out){
            if(copy_flow[out] > 0){
               if(edge_cost[out] < opt_cost){
                  opt_cost = edge_cost[out];
                  opt_arc = out;
               }
//               bottle_neck = std::max(bottle_neck, copy_flow[out]);
//               break;
            }
         }
         cur_node = g.target(opt_arc);
         path.push_back(opt_arc);
      }
      for(auto edge: path){
         --copy_flow[edge];
      }
      paths.push_back(path);
   }
}

void FlowNetwork::remove_low_cov_exon(const int cluster_left, const std::vector<float>& exon_doc,
                            std::list<std::pair<uint,uint>>& exon_boundaries)
{
  auto it = exon_boundaries.begin();
  while(it != exon_boundaries.end()){
     assert(it->second > it->first);
     auto it_start = next(exon_doc.begin(),it->first - cluster_left);
     auto it_end = next(exon_doc.begin(), it->second - cluster_left);
     double cov = accumulate(it_start, it_end, 0.0);
     cov = cov / (it->second - it->first);
     if(cov < kMinExonDoc)
       it =  exon_boundaries.erase(it);
     else if(it->second - it->first + 1 < kMinExonLen && cov < 5.0)
        it =  exon_boundaries.erase(it);
     else
        ++it;
  }
}


void FlowNetwork::filter_exon_segs(const std::vector<std::pair<uint,uint>>& paired_bars,
                         std::list<std::pair<uint,uint>>& exon_boundaries)
{
/*
  * filter exon segments if it does not have intron supporting
* */

//#ifdef DEBUG
//   for(auto ex = exon_boundaries.cbegin(); ex != exon_boundaries.cend(); ++ex){
//      std::cout<<"exon seg: "<<ex->first <<"-"<<ex->second<<std::endl;
//   }
//   for(auto in = paired_bars.cbegin(); in != paired_bars.cend(); ++in){
//      std::cout<<"intron: "<<in->first<<"-"<<in->second<<std::endl;
//   }
//#endif
   std::vector<size_t> dropoff;
   std::vector<std::pair<uint, uint>> left_coords;
   std::vector<std::pair<uint, uint>> right_coords;
   std::vector<std::pair<uint, uint>> e_boundaries;
   for(auto i: exon_boundaries){
      e_boundaries.push_back(i);
   }

   for(uint i = 0; i < paired_bars.size(); ++i){
      left_coords.push_back(std::pair<uint, uint>(paired_bars[i].first, i));
      right_coords.push_back(std::pair<uint, uint>(paired_bars[i].second, i));
   }
   sort(left_coords.begin(), left_coords.end(), comp_lt_first);
   sort(right_coords.begin(), right_coords.end(), comp_lt_first);

   for(size_t ex = 0 ; ex != e_boundaries.size(); ++ex){

      bool no_intron_on_right = false;
      auto l = lower_bound(left_coords.begin(), left_coords.end(), std::pair<uint, uint>(e_boundaries[ex].second+1,0), comp_lt_first);
      if(l != left_coords.end() && l->first == e_boundaries[ex].second+1){
         uint right = paired_bars[l->second].second;

         if(!binary_search(e_boundaries.begin(), e_boundaries.end(), std::pair<uint, uint>(right+1,0), comp_lt_first))
            no_intron_on_right = true;

      }
      else{
         no_intron_on_right = true;
      }

      bool no_intron_on_left = false;
      auto r = lower_bound(right_coords.begin(), right_coords.end(), std::pair<uint,uint>(e_boundaries[ex].first-1,0), comp_lt_first);
      if( r != right_coords.end() && r->first == e_boundaries[ex].first-1){
         uint left = paired_bars[r->second].first;
         if(!binary_search(e_boundaries.begin(), e_boundaries.end(), std::pair<uint, uint>(0,left-1), comp_lt_second)){
            no_intron_on_left = true;
         }
      }
      else{
         no_intron_on_left = true;
      }
      if(no_intron_on_left && no_intron_on_right){
         if(ex == 0){
            if( e_boundaries[ex].second +1 != e_boundaries[ex+1].first){
               dropoff.push_back(ex);
            }
         }
         else if(ex == e_boundaries.size() -1 ){
            if( e_boundaries[ex-1].second +1 != e_boundaries[ex].first)
               dropoff.push_back(ex);
         }
         else{
            if(e_boundaries[ex].second + 1 != e_boundaries[ex+1].first ||
               e_boundaries[ex].first -1 != e_boundaries[ex-1].second){
               dropoff.push_back(ex);
            }
         }
      }
   }

//#ifdef DEBUG
//   for(auto i: exon_boundaries)
//      cout<<"left: "<<i.first<<" right: "<<i.second<<endl;
//#endif

   std::vector<std::list<std::pair<uint,uint>>::iterator> drops;
   for(auto d: dropoff){
      drops.push_back(next(exon_boundaries.begin(),d));

   }

   for(auto&d: drops){
      exon_boundaries.erase(d);
   }

}

void FlowNetwork::filter_intron(const std::vector<GenomicFeature> &exons,
         std::map<std::pair<uint,uint>, IntronTable> &intron_counter)
{
   auto it = intron_counter.begin();
   while(it != intron_counter.end()){
      const IntronTable &intron = it->second;
      auto e1 = lower_bound(exons.begin(), exons.end(), intron.left-1, search_right);
      auto e2 = lower_bound(exons.begin(), exons.end(), intron.right+1, search_left);
      if(e1 == exons.end() || e2 == exons.end()) {
         it = intron_counter.erase(it);
      }
      else{
        if(e1->right() != intron.left-1 || e2->left() != intron.right + 1){
           it = intron_counter.erase(it);
        }
        else{
           ++it;
        }
      }
   }
}

void FlowNetwork::splicingGraph(const RefID & ref_id, const int &left, const std::vector<float> &exon_doc,
      std::map<std::pair<uint,uint>, IntronTable> &intron_counter,
      std::vector<GenomicFeature> &exons)
{
   /*
    * create non overlapping exon segments which will be
    * used as nodes in createNetwork();
    */


   std::vector<std::pair<uint,uint>> paired_bars; // two end intron boundaries
   std::vector<std::pair<uint,bool>> single_bars; // single intron boundaries

   for(auto i= intron_counter.cbegin(); i != intron_counter.cend(); ++i){
      paired_bars.push_back(std::pair<uint, uint> (i->first.first, i->first.second));
      single_bars.push_back(std::pair<uint, bool> (i->first.first, true));
      single_bars.push_back(std::pair<uint, bool> (i->first.second, false));
   }

   // unique element in single_bars
   sort(single_bars.begin(), single_bars.end(),comp_lt_first);
   auto newend = unique(single_bars.begin(), single_bars.end(),
         [](const std::pair<uint, bool> &lhs, const std::pair<uint, bool> &rhs){
         return (lhs.first == rhs.first && lhs.second == rhs.second) ;}
         );
   single_bars.erase(newend, single_bars.end());

   // unique element in paired_bars
   sort(paired_bars.begin(), paired_bars.end());
   auto new_end = unique(paired_bars.begin(), paired_bars.end());
   paired_bars.erase(new_end, paired_bars.end());
//#ifdef DEBUG
//   for(auto s: paired_bars)
//      cout<<"intorn: "<<s.first<<"-"<<s.second<<endl;
//   cout<<"---------------------"<<endl;
//#endif

   std::list<std::pair<uint,uint>> exon_boundaries;


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
         exon_boundaries.push_back(std::pair<uint,uint>(l,r));
         l = 0;
      }
   }
   if( l != 0 && l < left+exon_doc.size() )
      exon_boundaries.push_back(std::pair<uint,uint>(l, left+exon_doc.size()-1));


   /*
    * When some exonic coverage gaps exist due to low sequncing coverages.
    * This loop tries to fill in gaps and bring together separated exons.
    */
   auto iTer = exon_boundaries.begin();
   while(true){
      uint head = iTer->second;
      ++iTer;
      if(iTer == exon_boundaries.end()) break;
      uint tail = iTer->first;
      bool no_intron_overlap = true;
      bool no_intron_support = true;
      for(auto i= intron_counter.cbegin(); i != intron_counter.cend(); ++i){
         if(i->first.first <= tail && head <= i->first.second){
            no_intron_overlap = false;
         }
         if(i->first.first == head + 1 && tail-1 == i->first.second){
            no_intron_support =false;
         }
      }
      if(no_intron_overlap){
         if(tail - head < kMaxCoverGap1){
            iTer--;
            uint newStart = iTer->first;
            exon_boundaries.erase(iTer++);
            iTer->first = newStart;
         }
      }
      else{
         if(no_intron_support && tail - head < kMaxCoverGap2){
            iTer--;
            uint newStart = iTer->first;
            exon_boundaries.erase(iTer++);
            iTer->first = newStart;
         }
      }
   };


//      for(auto ex: exon_boundaries){
//         std::cout<<ex.first<<"-"<<ex.second<<std::endl;
//      }

   /*
    * for single exon genes
    * */

   if(paired_bars.size() == 0){
      uint l = exon_boundaries.front().first;
      uint r = exon_boundaries.back().second;
      //cout<<exon_boundaries.size()<<endl;
      exons.push_back(GenomicFeature(Match_t::S_MATCH, l, r-l+1));
      compute_exon_doc(left, exon_doc, exons);
      return;
   }

//#ifdef DEBUG
//   for(auto const & e: exon_boundaries){
//      std::cout<<"exon segments: "<<e.first<<"-"<<e.second<<std::endl;
//   }
//   for(auto const &s : single_bars){
//      std::cout<<s.first<<":"<<s.second<<"\t";
//   }
//   std::cout<<std::endl;
//#endif

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
            exon_boundaries.insert(++e, std::pair<uint,uint>(bar, temp));
         }
         else{
            e->second = bar;
            exon_boundaries.insert(++e, std::pair<uint,uint>(bar+1, temp));
         }
         --e;
         ++s;
      }
      else{
         ++e;
      }
   }

   /*
    * Filtering
    */
   auto it = exon_boundaries.begin();
   while(it != exon_boundaries.end()){
     if(it->second <= it->first){
        LOG_ERR("Unreal exon seg on ",ref_id,":", it->first, "-", it->second);
        //exit(0);
        it = exon_boundaries.erase(it);
     }
     else{
        ++it;
     }
   }
   //remove_low_cov_exon(left, exon_doc, exon_boundaries);
   filter_exon_segs(paired_bars, exon_boundaries);
   for(auto i: exon_boundaries){
      if(i.second - i.first +1 > 0)
         exons.push_back(GenomicFeature(Match_t::S_MATCH, i.first, i.second-i.first+1));
//      double sum = accumulate(exon_doc.begin() + i.first - left, exon_doc.begin() + i.second+1-left,0);
//      double avg_doc = sum/(i.second - i.first +1);
//      exons.back().avg_doc(avg_doc);
   }
   sort(exons.begin(), exons.end());
   compute_exon_doc(left, exon_doc, exons);
   filter_intron(exons, intron_counter);

}

bool FlowNetwork::createNetwork(
      const std::vector<Contig> &hits,
      const std::vector<GenomicFeature> &exons,
      const std::map<std::pair<uint,uint>, IntronTable> &intron_counter,
      const std::vector<std::vector<size_t>> &constraints,
      Graph::NodeMap<const GenomicFeature*> &node2feat,
      Graph::ArcMap<int> &cost_map,
      Graph::ArcMap<int> &min_flow_map,
      std::vector<std::vector<Graph::Arc>> &path_cstrs)
{
//#ifdef DEBUG
//   std::cout<<std::endl;
//   for(auto e: exons){
//      std::cout<<"exon: "<<e.left()<<"-"<<e.right()<<std::endl;
//   }
//   for(auto i:intron_counter){
//       const IntronTable &intron = i.second;
//       std::cout<<"intron: "<<intron.left<< "-" << intron.right<<std::endl;
//   }
//#endif
   assert(!hits.empty());
   if(exons.size() == 1) {
      //std::cout<<exons[0].left()<<"-"<<exons[0].right()<<std::endl;
      return true;
   }

   std::vector<Graph::Node> nodes;
   std::vector<Graph::Arc> arcs;
   std::map<const GenomicFeature*, Graph::Node> feat2node;
   for(size_t i = 0; i< exons.size(); ++i){
      nodes.push_back(_g.addNode());
      node2feat[nodes[i]] = &exons[i];
      feat2node[&exons[i]] = nodes[i];
   }
   if(exons.empty() || intron_counter.empty()) return false;
   // single exon case;

   // 1) add arc defined by intron
   for(auto i = intron_counter.cbegin(); i != intron_counter.cend(); ++i){
      const IntronTable &intron = i->second;
      auto e1 = lower_bound(exons.begin(), exons.end(), intron.left-1, search_right);
      auto e2 = lower_bound(exons.begin(), exons.end(), intron.right+1, search_left);
      if(e1 == exons.end() || e2 == exons.end()) {
#ifdef DEBUG
         //continue;
         assert(false);
         //std::cout<<intron.left<<"-"<<intron.right<<std::endl;
         //std::cout<<e1->left()<<"-"<<e1->right()<<std::endl;
#endif
         //return false;
      }

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

    /*
    * Single-end path constraint algorithm.
    */

   for(auto c: constraints){
      std::vector<Graph::Arc> path_cstr;

      size_t s_idx = c.front();
      size_t t_idx = c.back();
      const Graph::Node &s = feat2node[&exons[s_idx]];
      const Graph::Node &t = feat2node[&exons[t_idx]];

      //At least one of the internal nodes have >1 out-degree and >1 in-degree
      bool is_valid = false;
      for(size_t idx = 1; idx<c.size()-1; ++idx){
         const Graph::Node &n = feat2node[&exons[c[idx]]];
         if(inDeg[n] > 1 && outDeg[n] > 1)
            is_valid = true;
      }

      if(ArcLookUp<ListDigraph>(_g)(s,t) == INVALID && is_valid){ // no existing edge and valid constraint

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
         if(!path_cstr.empty())
            path_cstrs.push_back(path_cstr);
      }
   }


   /*
    * Paired-end path constraint algorithm. Currently defunct.
    */

//   for(auto c: constraints){
//      std::vector<Graph::Arc> path_cstr;
//
//      for(size_t i = 0; i< c.size()-1; ++i){
//         const Graph::Node &pre = feat2node[&exons[c[i]]];
//         const Graph::Node &sec = feat2node[&exons[c[i+1]]];
//
//         Graph::Arc arc_found = ArcLookUp<ListDigraph>(_g)(pre, sec);
//         if(arc_found == INVALID){
//            Dijkstra<Graph, Graph::ArcMap<int>> shortest_path(_g, cost_map);
//            shortest_path.run(pre);
//            std::cout << "The distance of node t from node s: "
//              << shortest_path.dist(sec) << std::endl;
//            if(shortest_path.dist(sec) == 0){
//               LOG("Path Constraints:\t");
//               LOG("(", exons[c[i]].left(), ",",exons[c[i]].right(),\
//             ")-(", exons[c[i+1]].left(),",",exons[c[i+1]].right(),")\t");
//            }
//            else{
//                std::stack<Graph::Arc> arcs_stack;
//                for (Graph::Node v=sec;v != pre; v=shortest_path.predNode(v)) {
//                   assert(v != INVALID);
//                   arcs_stack.push(shortest_path.predArc(v));
//                }
//                while(!arcs_stack.empty()){
//
//                   auto it = find(path_cstr.begin(), path_cstr.end(), arcs_stack.top());
//                   if(it != path_cstr.end());
//                      path_cstr.push_back(arcs_stack.top());
//                   arcs_stack.pop();
//                }
//            }
//
//         }
//         else{
//            auto it = find(path_cstr.begin(), path_cstr.end(), arc_found);
//            if( it != path_cstr.end());
//               path_cstr.push_back(arc_found);
//         }
//      }
//      if(path_cstr.size() > 1)
//         path_cstrs.push_back(path_cstr);
//   }


   if(path_cstrs.empty()){
      for(auto arc: arcs){
         min_flow_map[arc] = 1;
      }
      return true;
   }


   /*
    * Convert edges to path constraints while avoiding
    * duplications ( path constraints already contains edges)
    */
   std::vector<Graph::Arc> one_d_path_cstrs;
   for(auto p: path_cstrs){
     for(auto e: p)
        one_d_path_cstrs.push_back(e);
   }
   sort(one_d_path_cstrs.begin(), one_d_path_cstrs.end());
   auto new_end = unique(one_d_path_cstrs.begin(), one_d_path_cstrs.end());
   one_d_path_cstrs.erase(new_end, one_d_path_cstrs.end());

   for(auto edge: arcs){
      if(find(one_d_path_cstrs.begin(), one_d_path_cstrs.end(), edge) == one_d_path_cstrs.end()){
         std::vector<Graph::Arc> path;
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

void FlowNetwork::addWeight(const std::vector<Contig> &hits,
      const std::map<std::pair<uint,uint>, IntronTable> &intron_counter,
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
      if(arc_e - arc_s == 1){// if exon segs are next to each other
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
               num_read_support = i->second.total_junc_reads*kIntronEdgeWeight;
               break;
            }
         }
      }
      _max_weight = std::max(_max_weight, num_read_support);
      arc_map[arc] = num_read_support;
   }

   for(Graph::ArcIt arc(_g); arc != lemon::INVALID; ++arc){
      arc_map[arc] = _max_weight - arc_map[arc];
      assert(arc_map[arc] >= 0.0);
   }
}


   /*
    * Paired-end path constraint algorithm. Currently defunct.
    */
//std::vector<std::vector<size_t>> FlowNetwork::findConstraints(
//
//   const std::vector<GenomicFeature> &exons,
//   const std::vector<Contig> &hits)
//{
//   std::vector<std::vector<size_t>> result;
//   for(auto mp = hits.cbegin(); mp != hits.cend(); ++mp){
//      std::vector<size_t> constraint;
//      std::vector<GenomicFeature> introns;
//      for(auto feat = mp->_genomic_feats.cbegin(); feat != mp->_genomic_feats.cend(); ++feat){
//         if(feat->_match_op._code == Match_t::S_MATCH){
//            for(size_t i = 0; i< exons.size(); ++i){
//               if(GenomicFeature::overlap_in_genome(exons[i], feat->left(), feat->right() )  ){
//                  constraint.push_back(i);
//               }
//            }
//         }
//         // repeat Match_t::S_GAP
//         sort(constraint.begin(), constraint.end());
//         auto new_end = unique(constraint.begin(), constraint.end());
//         constraint.erase(new_end, constraint.end());
//         if(constraint.size() > 1){
//            result.push_back(constraint);
//         }
//      }
//   }
//   sort(result.begin(), result.end());
//   auto new_end = unique(result.begin(), result.end());
//   result.erase(new_end, result.end());
////   std::cout<<"num exons: "<<exons.size()<<std::endl;
////   std::cout<<"num constriants: "<<result.size()<<std::endl;
//   return result;
//}


   /*
    * Single-end path constraint algorithm.
    */

std::vector<std::vector<size_t>> FlowNetwork::findConstraints(

   const std::vector<GenomicFeature> &exons,
   const std::vector<Contig> &hits)
{
   std::vector<std::vector<size_t>> result;
   for(auto mp = hits.cbegin(); mp != hits.cend(); ++mp){
      std::vector<GenomicFeature> first_read;
      std::vector<GenomicFeature> second_read;

      bool is_first =true;
      for(auto feat = mp->_genomic_feats.cbegin(); feat != mp->_genomic_feats.cend(); ++feat){
         if(feat->_match_op._code == Match_t::S_GAP){
            is_first = false;
            continue;
         }
         if(is_first){
            first_read.push_back(*feat);
         }
         else{
            second_read.push_back(*feat);
         }
      }

      if(!first_read.empty()){
         std::vector<size_t> constraint;
         constraints_4_single_end(first_read, exons, constraint);
         if(constraint.size() > 2){
            result.push_back(constraint);
         }
      }
      if(!second_read.empty()){
         std::vector<size_t> constraint;
         constraints_4_single_end(first_read, exons, constraint);
         if(constraint.size() > 2){
            result.push_back(constraint);
         }
      }
   }
   sort(result.begin(), result.end());
   auto new_end = unique(result.begin(), result.end());
   result.erase(new_end, result.end());
   return result;
}


bool FlowNetwork::solveNetwork(const Graph::NodeMap<const GenomicFeature*> &node_map,
      const std::vector<GenomicFeature> &exons,
      const std::vector<std::vector<Graph::Arc>> &path_cstrs,
      Graph::ArcMap<int> &cost_map,
      Graph::ArcMap<int> &min_flow_map,
      std::vector<std::vector<GenomicFeature>> &transcripts)
{

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
   //if(exons[0].left() == 18415940){
#ifdef DEBUG
      for(Graph::NodeIt n(_g); n != lemon::INVALID; ++n){
         if( n == _source|| n == _sink) continue;
         std::cout<<_g.id(n)<<":"<<node_map[n]->left()<<"-"<<node_map[n]->right()<<std::endl;
      }
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
   std::vector<std::vector<Graph::Arc>> paths;
   flowDecompose(_g, flow, cost_map, _source, _sink, paths);

   //put paths into std::vector<std::vector<GenomicFeature>>
   //int i=0;
   for(auto p: paths){
      //cout<<"number "<<i<<endl;
      //++i;
      std::vector<GenomicFeature> tscp;
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
                  //cout<<" inside constraint exon: "<< node_map[n1]->left()<<"-"<<node_map[n1]->right()<<endl;
                  tscp.push_back(*node_map[n1]);
                  //if(idx+1 < cstr.size()-1)
                     if(node_map[n2]->left()-node_map[n1]->right() > 1){
                        //cout<<" inside constriant intron "<<node_map[n1]->right()+1<<"-"<<node_map[n2]->left()<<endl;
                        tscp.push_back(GenomicFeature(Match_t::S_INTRON, node_map[n1]->right()+1, node_map[n2]->left()-1-node_map[n1]->right()));
                     }
               }

               const Graph::Node &n1 = _g.source(cstr.back());
               const Graph::Node &n2 = _g.target(cstr.back());
               tscp.push_back(*node_map[n1]);
               if(node_map[n2]->left()-node_map[n1]->right() > 1){
                  tscp.push_back(GenomicFeature(Match_t::S_INTRON, node_map[n1]->right()+1, node_map[n2]->left()-1-node_map[n1]->right()));
               }

               break;
            }
         }
         // else it is edge originally
         if(is_edge){
            //cout<<node_map[s]->left()<<"-"<<node_map[s]->right()<<"\t";
            tscp.push_back(*node_map[arc_s]);
            if( i+1 < p.size()){
               if(node_map[arc_t]->left()-node_map[arc_s]->right() > 1){
                  tscp.push_back(GenomicFeature(Match_t::S_INTRON, node_map[arc_s]->right()+1, node_map[arc_t]->left()-1-node_map[arc_s]->right()));
               }
            }
         }
      }
      transcripts.push_back(tscp);
   }

   filter_short_transcripts(transcripts);
   if(transcripts.empty()) return false;

   return true;
}

void FlowNetwork::filter_short_transcripts(std::vector<std::vector<GenomicFeature>> &transcripts)
{
   for(auto it=transcripts.begin(); it != transcripts.end();){
      int len = 0;
      for(auto const& gf: *it){
         if(gf._match_op._code == Match_t::S_MATCH){
            len += gf._match_op._len;
         }
      }
      if(len < kMinTransLen ){
         it = transcripts.erase(it);
      }
      else{
         ++it;
      }
   }
}
