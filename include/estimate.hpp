/*
 * qp.h
 *
 *  Created on: Oct 15, 2015
 *      Author: ruolin
 */

#ifndef QP_H_
#define QP_H_
#include <contig.h>
#include <read.hpp>
class FaSeqGetter;
class Isoform{

public:
   Contig _contig;
   std::vector<GenomicFeature> _exon_segs;
   int _gene_id;
   int _isoform_id;
   std::string _isoform_str;
   std::string _gene_str;
   double _bais_factor;
   double _TPM = 0.0;
   double _FPKM = 0.0;
   std::string _TPM_s = "nan";
   std::string _FPKM_s = "nan";
   //Isoform() = default;
   Isoform(const std::vector<GenomicFeature>&, Contig, int, int);
   Isoform(const std::vector<GenomicFeature>&, Contig, std::string, std::string, int, int );
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
   std::set<std::pair<uint,uint>> _coords;
   double _GC_content;
public:
   std::set<Contig> _frags;
   void add_read_mass(float mass);

   std::map<int, double> _bin_weight_map; // iso -> bin_weight
   std::map<int, std::vector<std::pair<int,float>>> _iso_2_frag_lens; // iso -> (frag_len, mass)
   ExonBin(std::set<std::pair<uint,uint>> coordinates);
   uint left() const;
   uint right() const;
   float read_count() const;
   void add_frag_len(const int iso, const int frag_len, const float mass);
   double read_depth() const;
   int bin_len() const;
   std::vector<uint> bin_under_iso(const Isoform& iso,
         std::vector<std::pair<uint, uint>> & exon_coords) const;

   int effective_len(const std::vector<uint> & exons,
         const std::vector<uint>& implicit_idx,
         const int fl,
         const int rl
         ) const;
   bool operator==(const ExonBin& rhs) const;
   bool add_frag(const Contig& fg);
   int num_exons() const;
   int left_exon_len() const;
   double bin_gc_content(const std::shared_ptr<FaSeqGetter> &fa_getter, const int readlen) const;
   double bin_gc_content(const std::shared_ptr<FaSeqGetter> &fa_getter) const;
   double avg_frag_len() const;
   RefID ref_id() const;
};


class Estimation {
//typedef CGAL::Quadratic_program<double> Program;
//typedef CGAL::Quadratic_program_solution<ET> Solution;
   static const double _kMinTPM;
   static constexpr double _kMinFPKM = 1e-2;
   std::shared_ptr<InsertSize> _insert_size_dist;
   int _read_len;
   FILE* _p_log_file;
   void set_maps(const int& iso_id, const int& fg_len,
                const float& mass,
                const Contig& read,
                const std::set<std::pair<uint,uint>>& coords,
                std::map<std::set<std::pair<uint,uint>>, ExonBin> & exon_bin_map,
                std::map<int, std::set<std::set<std::pair<uint,uint>>>> &iso_2_bins_map);

public:
   Estimation(std::shared_ptr<InsertSize> insert_size, int read_len, FILE* tracker);
   void test();
   void overlap_exons(const std::vector<GenomicFeature>& exons,
                     const Contig& read,
                     std::set<std::pair<uint,uint>> &coords);

   void assign_exon_bin(
      const std::vector<Contig> &hits,
      const std::vector<Isoform> &transcripts,
      const std::vector<GenomicFeature> & exon_segs,
      std::map<std::set<std::pair<uint,uint>>, ExonBin> & exon_bin_map,
      std::map<int, std::set<std::set<std::pair<uint,uint>>>> &iso_2_bins_map);

   void theory_bin_weight(const std::map<int, std::set<std::set<std::pair<uint,uint>>>> &iso_2_bins_map,
                          const std::map<int,int> &iso_2_len_map,
                          const std::vector<Isoform>& isoforms,
                          std::map<std::set<std::pair<uint,uint>>, ExonBin> & exon_bin_map
                          );

   void empirical_bin_weight(const std::map<int, std::set<std::set<std::pair<uint,uint>>>> &iso_2_bins_map,
                          const std::map<int,int> &iso_2_len_map,
                          const int m,
                          std::map<std::set<std::pair<uint,uint>>, ExonBin> & exon_bin_map);

   void calculate_raw_iso_counts(const std::map<int, std::set<std::set<uint>>> &iso_2_bins_map,
          const std::map<std::set<uint>, ExonBin> & exon_bin_map);

   bool estimate_abundances(std::map<std::set<std::pair<uint,uint>>, ExonBin> & exon_bin_map,
                     const double mass,
                     std::map<int, int>& iso_2_len_map,
                     std::vector<Isoform>& isoforms,
                     bool with_bias_correction,
                     const std::shared_ptr<FaSeqGetter> & fa_getter);
   void calculate_bin_bias( std::map<std::set<std::pair<uint,uint>>, ExonBin> & exon_bin_map,
                            const std::shared_ptr<FaSeqGetter> &fa_getter,
                            std::vector<std::vector<double>> &bias);
};





class EmSolver{
   const double TOLERANCE = std::numeric_limits<double>::denorm_min();
   std::vector<double> _theta_after_zero;
   std::vector<int> _u; // observed data vector
   std::vector<std::vector<double>> _F; // sampling rate matrix, \alpha
   std::vector<std::vector<double>> _B; // bias matrix, \beta
   std::vector<std::vector<double>> _U; // hidden unobserved data matrix.
   static constexpr int _max_iter_num = 1000;
   static constexpr int _max_bias_it_num = 10;
   static constexpr int _max_theta_it_num = 5000;
   static constexpr int _max_out_it_num = 100;
   static constexpr double _theta_change_limit = 1e-2;
   static constexpr double _bias_change_limit = 1e-2;
public:
   std::vector<double> _theta;
   std::vector<double> _bias;
   EmSolver() = default;
   bool init( const int num_iso,
         const std::vector<int> &count,
         const std::vector<std::vector<double>> &model);

   bool init( const int num_iso,
         const std::vector<int> &count,
         const std::vector<std::vector<double>> &model,
         const std::vector<std::vector<double>> &bias);

   bool run();
};

int gap_ef(const int l_left, const int l_right, const int l_int, const int rl, const int gap);
int no_gap_ef(const int l_left, const int l_right, const int l_int, const int fl);
uint generate_pair_end(const Contig& ct, const Contig& s, int span,  SingleOrit_t orit, Contig & mp);

#endif /* QP_H_ */
