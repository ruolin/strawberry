#ifndef BIAS_HPP
#define BIAS_HPP
#include<vector>
#include "fasta.h"
#include "isoform.h"
#include "contig.h"
#include "alignments.h"

class LocusBias{
public:
   static std::vector<std::vector<double>> locus_bias(const Contig& , const std::vector<Isoform>& ) ; 
   static std::vector<double> iso_bias(const Contig&, const Isoform& ) ;

private:
   std::shared_ptr<HitCluster> _cluster;
   std::map<std::set<std::pair<uint, uint>>, ExonBin>& exon_bin_map;

};


#endif /* BIAS_HPP */
