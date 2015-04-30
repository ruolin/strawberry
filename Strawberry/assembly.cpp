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
//#ifdef DEBUG
//      cout<<"left: "<<l<<" right: "<<r<<endl;
//   cout<<"---------------------"<<endl;
//#endif
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

void FlowNetwork::createNetwork(
      const vector<Contig> &hits,
      const vector<GenomicFeature> &exons,
      const vector<IntronTable> &intron_counter,
      const vector<size_t> &bad_introns,
      const vector<vector<size_t>> &constraints,
      Graph::NodeMap<const GenomicFeature*> &node2feat,
      Graph::ArcMap<int> &cost_map,
      Graph::ArcMap<int> &min_flow_map,
      vector<vector<Graph::Arc>> &path_cstrs)
{

   vector<Graph::Node> nodes;
   vector<Graph::Arc> arcs;
   map<const GenomicFeature*, Graph::Node> feat2node;
   for(size_t i = 0; i< exons.size(); ++i){
      nodes.push_back(_g.addNode());
      node2feat[nodes[i]] = &exons[i];
      feat2node[&exons[i]] = nodes[i];
   }

   // 1) add arc defined by intron
   for(size_t i =0; i< intron_counter.size(); ++i){
      if( binary_search(bad_introns.begin(), bad_introns.end(), i) )
         continue;
      const IntronTable &intron = intron_counter[i];
      auto e1 = lower_bound(exons.begin(), exons.end(), intron.left-1, search_right);
      auto e2 = lower_bound(exons.begin(), exons.end(), intron.right+1, search_left);
      //cout<< intron.left << " vs " << intron.right<<endl;
      //cout<< e1->_genomic_offset+e1->_match_op._len<<" and "<< e2->_genomic_offset<<endl;
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
            assert(arc_found != INVALID);
            path_cstr.push_back(arc_found);
         }
         path_cstrs.push_back(path_cstr);
      }
   }

   if(path_cstrs.empty()){
      for(auto arc: arcs){
         min_flow_map[arc] = 1;
      }
      return;
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

void FlowNetwork::addWeight(const vector<Contig> &hits,
      const vector<IntronTable> &intron_counter,
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

vector<vector<size_t>> FlowNetwork::findConstraints(

   const std::vector<GenomicFeature> &exons,
   const std::vector<Contig> &hits)
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

void FlowNetwork::solveNetwork(const Graph::NodeMap<const GenomicFeature*> &node_map,
      const vector<GenomicFeature> &exons,
      const vector<vector<Graph::Arc>> &path_cstrs,
      Graph::ArcMap<int> &cost_map,
      Graph::ArcMap<int> &min_flow_map,
      vector<vector<GenomicFeature>> &transcripts){
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
      fprintf(stderr, "Infeasible or unbounded FlowNetwork flow\n");
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
}
