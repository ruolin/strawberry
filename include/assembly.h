/*
 * Builidng splicing graphs.
 * Solve the CMPC problem.
 */


#ifndef ASSEMBLY_H_
#define ASSEMBLY_H_
#include<vector>
#include <lemon/list_graph.h>
#include<map>
#include "common.h"
#include "utils.h"
#include "contig.h"


class IntronTable;
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
         const Graph::ArcMap<int> &cost,
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
//   void initGraph(const int &left,
//           const std::vector<float> &exon_doc,
//           const std::map<std::pair<uint,uint>, IntronTable> &intron_counter,
//           const std::vector<size_t> &bad_introns,
//           std::vector<GenomicFeature> &exons,
//           Graph::NodeMap<const GenomicFeature*> &node2feat);

   void splicingGraph(const RefID & ref_id, const int &left, const std::vector<float> &exon_doc,
         std::map<std::pair<uint,uint>, IntronTable> &intron_counter,
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

   void filter_short_transcripts(std::vector<std::vector<GenomicFeature>> &transcripts);
   void filter_exon_segs(const std::vector<std::pair<uint,uint>>& paired_bars,
                         std::list<std::pair<uint,uint>>& exon_boundaries);
   void remove_low_cov_exon(const int cluster_left, const std::vector<float>& exon_doc,
                            std::list<std::pair<uint,uint>>& exon_boundaries);
   void filter_intron(const std::vector<GenomicFeature> &exons,
         std::map<std::pair<uint,uint>, IntronTable> &intron_counter);
};

void compute_exon_doc(const int left, const std::vector<float>& exon_doc, std::vector<GenomicFeature>& exons);

inline std::vector<Contig> assemble_2_contigs(const std::vector<std::vector<GenomicFeature>> & assembled_feats,
                                       const RefID & ref_id,
                                       const Strand_t & strand)
{
    std::vector<Contig> results;
    for(auto const &feats: assembled_feats){
        std::vector<GenomicFeature> merged_feats;
        GenomicFeature::mergeFeatures(feats, merged_feats);
        Contig assembled_transcript(ref_id, 0, strand, 1.0, merged_feats, true);
        if(assembled_transcript.avg_doc() < kMinDepth4Contig) {
            continue;
        }
        PushAndReturnIdx<Contig>(assembled_transcript, results);
    }
    return results;
}

#endif /* ASSEMBLY_H_ */
