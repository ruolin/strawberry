/*
 * assembly.h
 *
 *  Created on: Mar 18, 2015
 *      Author: ruolin
 */

#ifndef ASSEMBLY_H_
#define ASSEMBLY_H_
#include<vector>
#include <lemon/list_graph.h>
#include<map>
#include "common.h"

//namespace lemon{
//   class ListDigraph;
//}
class IntronTable;
class GenomicFeature;
class Contig;
class ExonBin;
class Isoform;
typedef lemon::ListDigraph Graph;

class FlowNetwork{
   float _max_weight = 0.0;
   static bool hasFlow(const Graph &g, const Graph::ArcMap<int> & flow, const Graph::Node node){
      for(Graph::OutArcIt out(g, node); out != lemon::INVALID; ++out){
         if (flow[out] > 0)
            return true;
      }
      return false;
   }

   static void add_sink_source(Graph &g, Graph::Node &source, Graph::Node &sink);

   static bool comp_lt_first(const std::pair<uint, uint> & lhs, const std::pair<uint, uint> &rhs);
   static bool comp_lt_second(const std::pair<uint, uint> & lhs, const std::pair<uint, uint> &rhs);

   static bool search_left(const GenomicFeature &lhs, const uint rhs);

   static bool search_right(const GenomicFeature & lhs, const uint rhs);

   static void flowDecompose(const Graph &g,
         const Graph::ArcMap<int> &flow,
         const Graph::Node &source,
         const Graph::Node &sink,
         std::vector<std::vector<Graph::Arc>> &paths );

   static void flow2Transcript(const Graph &g,
         const std::vector<std::vector<Graph::Arc>> &paths,
         std::vector<std::vector<GenomicFeature>> &transcripts){

   }
public:
   Graph _g;
   Graph::Node _source;
   Graph::Node _sink;
   void initGraph(const int &left,
           const std::vector<float> &exon_doc,
           const std::map<std::pair<uint,uint>, IntronTable> &intron_counter,
           const std::vector<size_t> &bad_introns,
           std::vector<GenomicFeature> &exons,
           Graph::NodeMap<const GenomicFeature*> &node2feat);

   void splicingGraph(const int &left, const std::vector<float> &exon_doc,
         const std::map<std::pair<uint,uint>, IntronTable> &intron_counter,
         std::vector<GenomicFeature> &exons);

   bool createNetwork(
         const std::vector<Contig> &hits,
         const std::vector<GenomicFeature> &exons,
         const std::map<std::pair<uint,uint>, IntronTable> &intron_counter,
         const std::vector<std::vector<size_t>> &constraints,
         Graph::NodeMap<const GenomicFeature*> &node_map,
         Graph::ArcMap<int> &cost_map,
         Graph::ArcMap<int> &min_flow_map,
         std::vector<std::vector<Graph::Arc>> &path_cstrs);

   void addWeight(const std::vector<Contig> &hits,
         const std::map<std::pair<uint,uint>, IntronTable> &intron_counter,
         const Graph::NodeMap<const GenomicFeature*> &node_map,
         Graph::ArcMap<int> &arc_map);

   // return the positions of exons in
   std::vector<std::vector<size_t>> findConstraints(
         const std::vector<GenomicFeature> &exons,
         const std::vector<Contig> &hits);

   bool solveNetwork(const Graph::NodeMap<const GenomicFeature*> &node_map,
         const std::vector<GenomicFeature> &exons,
         const std::vector<std::vector<Graph::Arc>> &path_cstrs,
         Graph::ArcMap<int> &cost_map,
         Graph::ArcMap<int> &min_map,
         std::vector<std::vector<GenomicFeature>> &transcripts);
};

void assemble_2_contigs(const std::vector<std::vector<GenomicFeature>>& assembled_feats,
                        const RefID & ref_id,
                        const Strand_t & strand,
                        std::vector<Contig>& transcript);

#endif /* ASSEMBLY_H_ */
