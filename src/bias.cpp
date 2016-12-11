#include "bias.hpp"
#include "contig.h"

//vector<vector<double>> Bias::locus_bias(const Contig &reads, const vector<Isoform> &isoforms)
//{
//   vector<vector<double>> bx;
//      for(auto i=isoforms.cbegin(); i != isoforms.cend(); ++i){
//
//      }
//   return bx;
//}
//
//
//vector<double> Bias::iso_bias(const Contig& read, const Isoform& isoform) {
//   vector<double> bx(2,0.0);
//   uint s = Contig::read_start_from_iso(isoform._contig, read);
//   uint l = Contig::fragment_len(read, isoform._contig);
//   bx[0] = (double) s;
//   bx[1] = (double) l;
////   cout<<"Bias start"<<endl;
////   for(auto i:bx){
////      cout<<i<<"\t";
////   }
////   cout<<endl;
//   return bx;
//}

