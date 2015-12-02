/*
 * qp.h
 *
 *  Created on: Oct 15, 2015
 *      Author: ruolin
 */

#ifndef QP_H_
#define QP_H_
#include <contig.h>
#include <read.h>

class Isoform{

public:
   Contig _contig;
   int _gene_id;
   int _isoform_id;
   double _bais_factor;
   double _abundance = 0.0;
   //Isoform() = default;
   Isoform(Contig contig, int gene, int isoform);

};

class ExonBin{
   /*
    * ExonBin is a set a continuous exons defined by an overlapping fragment.
    */
private:
   //uint _left_most;
   //uint _right_most;
   //RefID _ref_id;

   float _whole_read_mass = 0;
   /*  exon boundaries left1,right1,left2,right2,... */
   set<uint> _coords;
public:
   set<const Contig*> _frags;
   void add_read_mass(float mass);
   /* trainscript_id -> pair<frag_len, read_mass>  */
   map<int, double> _iso_bias_map;
   map<int, vector<pair<int,float>>> _iso_2_frag_lens;
   ExonBin(set<uint> coordinates);
   uint left_most() const;
   uint right_most() const;
   float read_count() const;
   void add_frag_len(const int iso, const int frag_len, const float mass);
   double read_depth() const;
   int bin_len() const;
   bool operator==(const ExonBin& rhs) const;
   bool add_frag(const Contig* fg);
};


class Estimation {
//typedef CGAL::Quadratic_program<double> Program;
//typedef CGAL::Quadratic_program_solution<ET> Solution;

   shared_ptr<InsertSize> _insert_size_dist;
   int _read_len;

   void set_maps(const int& iso_id, const int& fg_len,
                const float& mass,
                const Contig* read,
                const set<uint>& coords,
                map<set<uint>, ExonBin> & exon_bin_map,
                map<int, set<set<uint>>> &iso_2_bins_map);

public:
   Estimation(shared_ptr<InsertSize> insert_size, int read_len);
   void test();
   void overlap_exons(const vector<GenomicFeature> exons,
                     const uint left,
                     const uint right,
                     set<uint> &coords);

   void assign_exon_bin(
      const vector<Contig> &hits,
      const vector<Isoform> &transcripts,
      const vector<GenomicFeature> & exon_segs,
      map<set<uint>, ExonBin> & exon_bin_map,
      map<int, set<set<uint>>> &iso_2_bins_map);

   void calcuate_bin_bias(const map<int, set<set<uint>>> &iso_2_bins_map,
                          const map<int,int> &iso_2_len_map,
                          const int m,
                          map<set<uint>, ExonBin> & exon_bin_map);

   void calculate_raw_iso_counts(const map<int, set<set<uint>>> &iso_2_bins_map,
          const map<set<uint>, ExonBin> & exon_bin_map);

   void estimate_abundances(const map<set<uint>, ExonBin> & exon_bin_map,
                     const double mass,
                     vector<Isoform>& isoforms,
                     bool use_qp,
                     bool with_bias_correction = true);
};


uint generate_pair_end(const Contig& ct, const uint& start, int span,  SingleOrit_t orit);







#endif /* QP_H_ */
