#ifndef BIAS_HPP
#define BIAS_HPP
#include<vector>


class Isoform;
class Contig;
class Bias{
public:
   static std::vector<std::vector<double>> locus_bias(const Contig& , const std::vector<Isoform>& ) ; 
   static std::vector<double> iso_bias(const Contig&, const Isoform& ) ; 
};


#endif /* BIAS_HPP */
