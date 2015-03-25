/*
 * assembly.h
 *
 *  Created on: Mar 18, 2015
 *      Author: ruolin
 */

#ifndef ASSEMBLY_H_
#define ASSEMBLY_H_
#include<vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include<boost/version.hpp>
#if (BOOST_VERSION < 103800)
#include <boost/vector_property_map.hpp>
#else
#include<boost/property_map/vector_property_map.hpp>
#endif

typedef boost::adjacency_list<boost::vecS,
                             boost::vecS,
                             boost::bidirectionalS,
                             boost::property<boost::vertex_distance_t, int>,
                             boost::property<boost::edge_weight_t, int>> DAG;

typedef boost::graph_traits<DAG>::vertex_descriptor DAGVertex;
//typedef boost::property_map<DAG, boost::vertex_name_t>::type name_map_t;
//typedef boost::property_map<DAG, boost::edge_weight_t>::type edge_weight_map_t;

class IntronTable;
class GenomicFeature;

void assemble(const int &left,
        const std::vector<float> &exon_doc,
        const std::vector<IntronTable> &intron_counter,
        const std::vector<size_t> &bad_introns,
        std::vector<GenomicFeature> &exons,
        std::vector<GenomicFeature> &introns);

void exonDoc2GFeats(const int &left, const std::vector<float> &exon_doc,
      const std::vector<IntronTable> &intron_counter, const std::vector<size_t> &bad_introns,
      std::vector<GenomicFeature> &exons);

void flowCycleCanceling();
bool createDAG();
#endif /* ASSEMBLY_H_ */
