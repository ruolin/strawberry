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
#include <boost/graph/dag_shortest_paths.hpp>
#include <iostream>

#include <lemon/smart_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>

//#include <lemon/preflow.h>

using namespace lemon;
typedef ListDigraph Graph;
typedef int LimitValueType;


void assemble( const int &left,
        const std::vector<float> &exon_doc,
        const std::vector<IntronTable> &intron_counter,
        const std::vector<size_t> &bad_introns,
        std::vector<GenomicFeature> &exons,
        std::vector<GenomicFeature> &introns)
{
   exonDoc2GFeats(left, exon_doc, intron_counter, bad_introns, exons);
}

// comparator used for searching in the split_bars object, which records the
// intron boundaries. We use pair<uint,bool> to indicate its position and
// whether it is left boundary or right boundary.
// left = true; right = false

bool comp_lt(const pair<uint, bool> & lhs, const pair<uint,bool> &rhs){
   return lhs.first < rhs.first;
}

void exonDoc2GFeats(const int &left, const std::vector<float> &exon_doc,
      const std::vector<IntronTable> &intron_counter, const std::vector<size_t> &bad_introns,
      std::vector<GenomicFeature> &exons)
{

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
               e->second == next(e)->first-1) ;
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
         if( ( binary_search(split_bars.begin(), split_bars.end(), pair<uint, bool>(e->first-1,true), comp_lt) ||
               e->first == prev(e)->second+1
             ) &&
             (
               binary_search(split_bars.begin(), split_bars.end(), pair<uint, bool>(e->second+1,true), comp_lt) ||
               e->second == next(e)->first-1
             ));
         else
            dropoff.push_back(e);
      }
   }
   for(auto d: dropoff)
      exon_boundries.erase(d);
#ifdef DEBUG
   for(auto i: exon_boundries){
      cout<<"left: "<<i.first<<" right: "<<i.second<<endl;
   }
#endif
}






void flowCycleCanceling()
{
   /**
    * calculate min cost max flow via network canceling
    * sample graph (arcs are dircted from left to right
    *
    *     2     5
    *    / \   / \
    *   /   \ /   \
    *  0     4--6--1
    *   \   / \   /
    *    \ /   \ /
    *     3     7
    *
    * vertex 0 represents the source, 1 target
    * cost are set to 1, ecept for source and target arcs
    * cost for (3,4) = 2
    */

   Graph g;

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g), c_i(g);

   // source and sink vertices
   Graph::Node n0 = g.addNode(), n1 = g.addNode();

   // other nodes/vertices
   Graph::Node n2 = g.addNode(), n3 = g.addNode(), n4 = g.addNode(), n5 = g.addNode(),
      n6 = g.addNode(), n7 = g.addNode();

   Graph::Arc a02 = g.addArc(n0, n2), a03 = g.addArc(n0, n3), a24 = g.addArc(n2, n4),
         a34 = g.addArc(n3, n4), a45 = g.addArc(n4, n5), a46 = g.addArc(n4, n6),
         a47 = g.addArc(n4, n7), a51 = g.addArc(n5, n1), a61 = g.addArc(n6, n1),
         a71 = g.addArc(n7, n1);

   std::vector<Graph::Arc> arcs = { a02, a03, a24, a34, a45, a46, a47, a51, a61, a71 };

        // arcs for circulation! (source to target and vice versa)
   Graph::Arc ts = g.addArc(n1, n0);
   Graph::Arc st = g.addArc(n0, n1);


   // graph must be finished before initializing network canceling !
   NetworkSimplex<Graph, LimitValueType, LimitValueType> network(g);

        // bounds for circulation arcs
   l_i[ts] = 0; u_i[ts] = network.INF; c_i[ts] = 0;
   l_i[st] = 0; u_i[ts] = network.INF; c_i[ts] = 0;

        // set cost for arcs

   for ( auto & a : arcs )
   {
      l_i[a] = 0;
      u_i[a] = network.INF;
      c_i[a] = (g.source(a) == n0|| g.target(a) == n1 ) ? 0 : 1;
   }

   for (int i=4; i<8; ++i){
      l_i[arcs[i]] = 1;
   }
   c_i[a34] = 2; // arc (3,4)

        //set lower/upper bounds, cost
   network.lowerMap(l_i).upperMap(u_i).costMap(c_i);

   NetworkSimplex<Graph>::ProblemType ret = network.run();

   // get flow of arcs
   ListDigraph::ArcMap<LimitValueType> flow(g);
   network.flowMap(flow);

   digraphWriter(g).                  // write g to the standard output
           arcMap("cost", c_i).          // write 'cost' for for arcs
           arcMap("flow", flow).          // write 'flow' for for arcs
           arcMap("l_i", l_i).
           arcMap("u_i", u_i).
           node("source", n0).             // write s to 'source'
           node("target", n1).             // write t to 'target'
           run();

   switch ( ret )
   {
   case NetworkSimplex<Graph>::INFEASIBLE:
      std::cerr << "INFEASIBLE" << std::endl;
   break;
   case NetworkSimplex<Graph>::OPTIMAL:
      std::cerr << "OPTIMAL" << std::endl;
   break;
   case NetworkSimplex<Graph>::UNBOUNDED:
      std::cerr << "UNBOUNDED" << std::endl;
   break;
   }

   // determine flow from source and flow to target
   double flowToT = 0, flowToS = 0;
   for (ListDigraph::ArcMap<LimitValueType>::ItemIt a(flow); a != INVALID; ++a)
   {
      if ( g.target(a) == n1 )
         flowToT += network.flow(a);

      if ( g.source(a) == n0 )
         flowToS += network.flow(a);
   }

   std::cerr << "flow to t = " << flowToT << " flow from s = " << flowToS << std::endl;
   std::cerr << "total cost = " << network.totalCost<double>() << std::endl;

}


