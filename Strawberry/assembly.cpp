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
#include <lemon/bfs.h>

using namespace lemon;
typedef int LimitValueType;

bool FlowNetwork::comp_lt(const std::pair<uint, bool> & lhs, const std::pair<uint,bool> &rhs){
      return lhs.first < rhs.first;
   }

bool FlowNetwork::search_left(const GenomicFeature &lhs, const uint rhs){
   return lhs._genomic_offset < rhs;
}

bool FlowNetwork::search_right(const GenomicFeature & lhs, const uint rhs){
   uint right = lhs._genomic_offset + lhs._match_op._len -1;
   return right < rhs;
}

void FlowNetwork::initGraph(const int &left,
        const vector<float> &exon_doc,
        const vector<IntronTable> &intron_counter,
        const vector<size_t> &bad_introns,
        vector<GenomicFeature> &exons,
        Graph::NodeMap<const GenomicFeature*> &node2feat)
{
   splicingGraph(left, exon_doc, intron_counter, bad_introns, exons);
   createNetwork(exons, intron_counter, bad_introns, node2feat);
}

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
      const std::vector<IntronTable> &intron_counter, const std::vector<size_t> &bad_introns,
      std::vector<GenomicFeature> &exons)
{
   /*
    * create non overlapping exon segments which will be
    * used as nodes in createNetwork();
    */
   vector<pair<uint,bool>> split_bars;
   for(size_t i=0; i<intron_counter.size(); ++i){
      if(binary_search(bad_introns.begin(), bad_introns.end(), i)){
         LOG("The filtered intron ", intron_counter[i].left, "-", intron_counter[i].right, " is not further used in building splicing graph");
         continue;
      }
      split_bars.push_back(pair<uint, bool> (intron_counter[i].left, true));
      split_bars.push_back(pair<uint, bool> (intron_counter[i].right, false));
   }
   sort(split_bars.begin(), split_bars.end(),comp_lt);
   auto new_end = unique(split_bars.begin(), split_bars.end(),
         [](const pair<uint, bool> &lhs, const pair<uint, bool> &rhs){
         return lhs.first == rhs.first;}
         );
   split_bars.erase(new_end, split_bars.end());
   list<pair<uint,uint>> exon_boundries;

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
         exon_boundries.push_back(pair<uint,uint>(l,r));
         l = 0;
      }
   }
   if( l != 0 && l < left+exon_doc.size() )
      exon_boundries.push_back(pair<uint,uint>(l, left+exon_doc.size()-1));

   /*
    * for single exon genes
    * */
   if(split_bars.size() == 0){
      uint l = exon_boundries.front().first;
      uint r = exon_boundries.back().second;
      exons.push_back(GenomicFeature(Match_t::S_MATCH, l, r-l+1));
#ifdef DEBUG
      cout<<"left: "<<l<<" right: "<<r<<endl;
   cout<<"---------------------"<<endl;
#endif
      return;
   }

   /*
    * further divided preliminary exon segments into smaller pieces based on intorn boundries.
    * */
   auto e = exon_boundries.begin();
   size_t s = 0;
   while(e != exon_boundries.end() && s < split_bars.size()){

      uint bar = split_bars[s].first;
      bool left = split_bars[s].second;
      if(bar < e->first){
         ++s;
      }
      else if(bar >= e->first && bar <= e->second){
         uint temp = e->second;
         if(left){
            e->second = bar-1;
            exon_boundries.insert(++e, pair<uint,uint>(bar, temp));
         }
         else{
            e->second = bar;
            exon_boundries.insert(++e, pair<uint,uint>(bar+1, temp));
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
   vector<list<pair<uint,uint>>::iterator> dropoff;
   for(e = exon_boundries.begin(); e!= exon_boundries.end(); ++e){
      if(*e == exon_boundries.front()){
         if(binary_search(split_bars.begin(), split_bars.end(), pair<uint, bool>(e->second+1,true), comp_lt) ||
               e->second == next(e)->first-1);
         else
            dropoff.push_back(e);

      }
      else if(*e == exon_boundries.back()){
         if( binary_search(split_bars.begin(), split_bars.end(), pair<uint, bool>(e->first-1,true), comp_lt) ||
               e->first == prev(e)->second+1);
         else
            dropoff.push_back(e);
      }
      else{

         // first two conditiosn cover situations when exon
         // segments are created by coverage gaps not introns.
         // The gap situation can be found in AT1G01020.
         // It can be addressed here. But we leave this problem to buildDAG();
         if( next(e)->first - e->second > 1 &&
             !binary_search(split_bars.begin(), split_bars.end(), pair<uint, bool>(e->second+1,true), comp_lt)
           )
         {
            if( binary_search(split_bars.begin(), split_bars.end(), pair<uint, bool>(e->first-1,true), comp_lt));
            else{
               dropoff.push_back(e);
            }
            continue;
         }

         if (e->first - prev(e)->second > 1 &&
               !binary_search(split_bars.begin(), split_bars.end(), pair<uint, bool>(e->first-1,true), comp_lt)
            )
         {
            if(binary_search(split_bars.begin(), split_bars.end(), pair<uint, bool>(e->second+1,true), comp_lt)) ;
            else{
               dropoff.push_back(e);
            }
            continue;
         }

         if( ( binary_search(split_bars.begin(), split_bars.end(), pair<uint, bool>(e->first-1,true), comp_lt) ||
               e->first == prev(e)->second+1
             ) &&
             (
               binary_search(split_bars.begin(), split_bars.end(), pair<uint, bool>(e->second+1,true), comp_lt) ||
               e->second == next(e)->first-1
             ));
         else{
            dropoff.push_back(e);
         }
      }

   }
   //cout<<"ha"<<exon_boundries.front().first<<" "<<exon_boundries.front().second<<endl;
   for(auto d: dropoff){
      exon_boundries.erase(d);
   }
   for(auto i: exon_boundries){
      exons.push_back(GenomicFeature(Match_t::S_MATCH, i.first, i.second-i.first+1));
   }
   sort(exons.begin(), exons.end());
#ifdef DEBUG
   for(auto i: exon_boundries){
      cout<<"left: "<<i.first<<" right: "<<i.second<<endl;
   }
   cout<<"---------------------"<<endl;
#endif
}

void FlowNetwork::createNetwork(const vector<GenomicFeature> &exons,
      const vector<IntronTable> &intron_counter,
      const vector<size_t> &bad_introns,
      Graph::NodeMap<const GenomicFeature*> &node2feat){

   vector<Graph::Node> nodes;
   //vector<Graph::Arc> arcs;
   map<const GenomicFeature*, Graph::Node> feat2node;
   for(size_t i = 0; i< exons.size(); ++i){
      nodes.push_back(_g.addNode());
      node2feat[nodes[i]] = &exons[i];
      feat2node[&exons[i]] = nodes[i];
   }

   for(size_t i =0; i< intron_counter.size(); ++i){
      if( binary_search(bad_introns.begin(), bad_introns.end(), i) )
         continue;
      const IntronTable &intron = intron_counter[i];
      auto e1 = lower_bound(exons.begin(), exons.end(), intron.left-1, search_right);
      auto e2 = lower_bound(exons.begin(), exons.end(), intron.right+1, search_left);
      //cout<< intron.left << " vs " << intron.right<<endl;
      //cout<< e1->_genomic_offset+e1->_match_op._len<<" and "<< e2->_genomic_offset<<endl;
      Graph::Arc a = _g.addArc( feat2node[&(*e1)],  feat2node[&(*e2)]);
      //arcs.push_back(a);
   }

   for(size_t i=0; i < exons.size()-1; ++i){
      uint right = exons[i]._genomic_offset + exons[i]._match_op._len -1;
      if(exons[i+1]._genomic_offset - right == 1){
         const Graph::Node &node_a = feat2node[&exons[i]];
         const Graph::Node &node_b = feat2node[&exons[i+1]];
         Graph::Arc a = _g.addArc(node_a, node_b);
         //arcs.push_back(a);
      }
   }
//   for(auto a: arcs){
//      l[a] = 1;
//   }

//    digraphWriter(g).                  // write g to the standard output
//           arcMap("cost", c).          // write 'cost' for for arcs
//           arcMap("flow", flow).          // write 'flow' for for arcs
//           arcMap("l_i", l).
//           arcMap("u_i", u).
//           node("source", source).             // write s to 'source'
//           node("target", sink).             // write t to 'target'
//           run();
//   for(Graph::ArcIt a(g); a != INVALID; ++a){
//      if(g.source(a) == source ){
//         //cout<<"haha"<<endl;
//         cout<< node2feat[g.target(a)]->_genomic_offset<<endl;
//      }
//      if(g.target(a) == sink){
//         cout<< node2feat[g.source(a)]->_genomic_offset<<endl;
//      }
//      cout<<l[a]<<endl;
//   }
}

void FlowNetwork::addWeight(const std::vector<Contig> &hits,
      const std::vector<IntronTable> &intron_counter,
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
         for(auto i : intron_counter){
            if(arc_s == i.left && arc_e == i.right){
               num_read_support = i.total_junc_reads;
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

void FlowNetwork::solveNetwork(Graph::ArcMap<int> &cost_map){
   add_sink_source(_g, _source, _sink);
   Graph::ArcMap<LimitValueType> l(_g), u(_g);
   NetworkSimplex<Graph, LimitValueType, LimitValueType> FlowNetwork(_g);

   for(Graph::ArcIt arc(_g); arc != INVALID; ++arc){
      u[arc] = FlowNetwork.INF;
      l[arc] = 1;
   }
   FlowNetwork.lowerMap(l).upperMap(u).costMap(cost_map);

   NetworkSimplex<Graph>::ProblemType ret = FlowNetwork.run();

   // get flow of arcs
   Graph::ArcMap<LimitValueType> flow(_g);
   FlowNetwork.flowMap(flow);
   if(ret == NetworkSimplex<Graph>::INFEASIBLE || ret == NetworkSimplex<Graph>::UNBOUNDED){
      fprintf(stderr, "Infeasible or unbounded FlowNetwork flow\n");
   }
   vector<vector<Graph::Arc>> paths;
   vector<vector<GenomicFeature>> transcripts;
   flowDecompose(_g, flow, _source, _sink, paths);
   for(auto p: paths){
      for(auto e: p){
         cout<<cost_map[e]<<"\t";
      }
      cout<<endl;
   }
}
