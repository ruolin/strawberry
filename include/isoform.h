//
// Created by ruolin on 12/10/16.
//

#ifndef STRAWBERRY_ISOFORM_H
#define STRAWBERRY_ISOFORM_H

#include <string>
#include <algorithm>
#include "contig.h"



class Isoform{

   int _isoform_id;
   int _gene_id;
public:
   Contig _contig;
   std::vector<GenomicFeature> _exon_segs; // non-overlapping exon segments
   std::string _isoform_str;
   std::string _gene_str;
   double _bais_factor;
   double _TPM = 0.0;
   double _FPKM = 0.0;
   std::string _TPM_s = "nan";
   std::string _FPKM_s = "nan";
   //Isoform() = default;
   Isoform(const std::vector<GenomicFeature>& exons, Contig contig,
           std::string gene_name, std::string iso_name, int gene_id):
           _contig(contig), _gene_str(gene_name), _isoform_str(iso_name),
           _gene_id(gene_id)
   {
      for(uint i = 0; i< exons.size(); ++i){
         if(Contig::is_compatible(_contig, exons[i])){
            _exon_segs.push_back(exons[i]);
         }
      }
      _bais_factor = 0.0;
   }

   Isoform(const std::vector<GenomicFeature>& exons, Contig contig, int gene):
           Isoform(exons, contig, "default gene", "default iso", gene)
   {
      for(uint i = 0; i< exons.size(); ++i){
         if(Contig::is_compatible(_contig, exons[i])){
            _exon_segs.push_back(exons[i]);
         }
      }
      _bais_factor = 0.0;
   }

   decltype(auto) id() { return (_isoform_id);}
   decltype(auto) id() const { return (_isoform_id);}
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
   int _ebidx;
   double _GC_content;

   int no_gap_ef(const int l_left, const int l_right, const int l_int, const int fl) const
/*
 * Effective length of a fragment hitting both ends but
 * do not consider whether it hits the inner segments.
 */
   {
      if(fl < l_int + 2) return 0;
      if(fl > l_left + l_right + l_int) return 0;
      int mid = fl - l_int -1;
      return std::min(l_left, mid) + std::min(l_right, mid) - mid;
   }

   int gap_ef(const int l_left, const int l_right, const int l_int, const int rl, const int gap) const
/*
 * Effective length of a fragment hitting both ends but
 * CANNOT hit any inner segmetns!
 * This time we need to consider the read gap between two ends.
 */
   {
      if(2*rl+gap < l_int +2) return 0;
      if(2*rl+gap > l_left + l_right + l_int) return 0;
      int start = std::max(rl, l_left + l_int - gap -1);
      int end = std::min(l_left, l_left+l_right+l_int - gap - rl);
      return ( std::max(0, end-start));
   }


public:
   std::set<Contig> _frags;
   void add_read_mass(float mass);

   std::map<int, double> _bin_weight_map; // iso -> bin_weight
   std::map<int, std::vector<std::pair<int,float>>> _iso_2_frag_lens; // iso -> (frag_len, mass)
   ExonBin(std::set<std::pair<uint,uint>> coordinates): _coords(coordinates){};
   uint left() const;
   uint right() const;
   float read_count() const;
   void add_frag_len(const int iso, const int frag_len, const float mass);
   double read_depth() const;
   int bin_len() const;

   decltype(auto) id() const {return (_ebidx);}
   decltype(auto) id() {return (_ebidx);}

   std::vector<uint> bin_under_iso(const Isoform& iso,
                                   std::vector<std::pair<uint, uint>> & exon_coords) const;

   bool is_compatible(const Isoform& iso) const {
      std::set<std::pair<uint,uint>> iso_exon_segs;
      for (const auto& ex : iso._exon_segs) {
         iso_exon_segs.emplace(std::make_pair<uint, uint>(ex.left(), ex.right()));
      }
      return std::includes(iso_exon_segs.begin(), iso_exon_segs.end(), _coords.begin(), _coords.end() );
   }

   int len_under_iso(const Isoform& iso) const {
      int result = 0;
      if (!is_compatible(iso)) return result;

      std::vector<std::pair<uint, uint>> coords;
      std::vector<uint> implicit_idx =bin_under_iso(iso, coords);
      if (!implicit_idx.empty()) return result;
      for (const auto& c: coords) result += c.second - c.first + 1;
      return result;
   }


   int effective_len(const std::vector<uint> & exons,
                     const std::vector<uint>& implicit_idx,
                     const int fl,
                     const int rl) const;

   bool add_frag(const Contig& fg);
   int num_exons() const;
   int left_exon_len() const;
   double bin_gc_content(const std::shared_ptr<FaSeqGetter> &fa_getter, const int readlen) const;
   double bin_gc_content(const std::shared_ptr<FaSeqGetter> &fa_getter) const;
   double avg_frag_len() const;
   RefID ref_id() const;

   friend bool operator==(const ExonBin& lhs, const ExonBin& rhs);
   friend bool operator!=(const ExonBin& lhs, const ExonBin& rhs);
   friend bool operator<(const ExonBin& lhs, const ExonBin& rhs);
   friend std::ostream& operator<<(std::ostream&, const ExonBin& );
};

   inline std::ostream& operator<<(std::ostream& os, const ExonBin& eb){
      for (auto c: eb._coords) {
         os <<"[" << c.first <<"-" <<c.second<<"]("<<c.second-c.first+1<<")\t";
      }
      os<<std::endl;
      return os;
   }

inline bool operator==(const ExonBin& lhs, const ExonBin& rhs) {
   return lhs._coords == rhs._coords;
}

inline bool operator!=(const ExonBin& lhs, const ExonBin& rhs) {
   return !(lhs == rhs);
}

inline bool operator< (const ExonBin& lhs, const ExonBin& rhs) {
   return lhs._coords < rhs._coords;
}


inline uint ExonBin::left() const
{
   return _coords.cbegin()->first;
}


inline uint ExonBin::right() const
{
   return _coords.crbegin()->second;
}


inline void ExonBin::add_read_mass(float mass)
{
   _whole_read_mass += mass;
}

inline bool ExonBin::add_frag(const Contig& fg)
{
   auto ret = _frags.insert(fg);
//#ifdef DEBUG
//   if(ret.second == false){
//      cout<<(*ret.first).left()<<" has existed"<<endl;
//
//   }
//   if(ret.second == true){
//      cout<<fg.left()<<" has been added"<<endl;
//   }
//#endif
   return true;
}

inline int ExonBin::num_exons() const
{
   return _coords.size();
}

inline float ExonBin::read_count() const
{
   float sum = 0.0 ;
   for(auto it = _frags.cbegin(); it != _frags.cend(); ++it ){
//#ifdef DEBUG
//      cout<<"it->left "<<it->left()<<" mass: "<<it->mass()<<endl;
//#endif
      sum += it->mass();
   }
   return sum;
}

inline int ExonBin::bin_len() const
{
   int bin_len = 0;
   for(auto c = _coords.cbegin(); c != _coords.cend(); ++c){
      bin_len += c->second - c->first +1;
   }
   return bin_len;
}

inline double ExonBin::bin_gc_content(const std::shared_ptr<FaSeqGetter> &fa_getter, const int readlen) const
{
   int start = 0;
   int len = 0;
   if(num_exons() == 1){
      start = _coords.cbegin()->first;
      len = bin_len();
      return gc_content( fa_getter->fetchSeq(start, len));

   }
   else if(num_exons() == 2){
      double total_gc = 0.0;
      for(auto it = _coords.cbegin(); it != _coords.cend(); ++it){
         start = std::max(it->first, it->second - readlen + 1);
         len = it->second  - start +1;
         total_gc += gc_content( fa_getter->fetchSeq(start, len));
      }
      return total_gc/2;
   }
   else{
      std::vector<double> gc_ratios;
      for(auto it = _coords.cbegin(); it != _coords.cend(); ++it){
         if( *it == *_coords.cbegin() || *it == *_coords.crbegin()) continue;
         start = it->first;
         len += it->second - it->first +1;
         gc_ratios.push_back( gc_content( fa_getter->fetchSeq(start, len))  );
      }
      auto it = std::max_element(gc_ratios.begin(), gc_ratios.end());
      return *it;
   }

}


inline double ExonBin::bin_gc_content(const std::shared_ptr<FaSeqGetter> &fa_getter) const
{
   std::vector<double> gc_ratios;
   for(auto it= _coords.cbegin(); it != _coords.cend(); ++it){
      double gc = gc_content(fa_getter->fetchSeq(it->first, it->second - it->first+1) );
      gc_ratios.push_back(gc);
   }
   double total_gc = std::accumulate(gc_ratios.begin(), gc_ratios.end(), 0.0);
   return total_gc/gc_ratios.size();
}

inline double ExonBin::avg_frag_len() const
{
   assert(_iso_2_frag_lens.empty() == false);
   float total_mass = 0;
   double total_len = 0;
   for(auto const &by_iso : _iso_2_frag_lens){
      for(auto const &f: by_iso.second){
         total_len += f.first * f.second;
         total_mass += f.second;
      }
   }
   return total_len / total_mass;
}

inline RefID ExonBin::ref_id() const
{
   assert(_frags.empty() == false);
   return _frags.cbegin()->ref_id();
}


inline std::vector<uint> ExonBin::bin_under_iso(const Isoform& iso,
                                    std::vector<std::pair<uint, uint>> & exon_coords) const
/*
 * NOTE: The number of segments of a bin might be less than the number of
 * segments under its compatible isoform. This is because the fragment gap.
 *
 * Return two things.
 * First, the implicit whole exon segments under a isoform. This is returned by second input argument.
 * Second, the idx of the exon segments that are implicit. This is returned by return value.
 */
{
   std::vector<uint> idx;
   const std::vector<GenomicFeature> & exons = iso._exon_segs;
   std::vector<uint> start_pos(exons.size());
   for(uint i=0; i<exons.size(); ++i){
      start_pos[i] = exons[i].left();
   }

   std::vector<uint>::const_iterator low = lower_bound(start_pos.cbegin(), start_pos.cend(), this->left());

   assert(low != start_pos.cend());
   std::vector<uint>::const_iterator up  = lower_bound(start_pos.cbegin(), start_pos.cend(), _coords.crbegin()->first );
   assert(up  != start_pos.cend());
   for(std::vector<uint>::const_iterator it = low; it != up; ++it){
      const GenomicFeature & cur_exon = exons[distance(start_pos.cbegin(),it)];
      exon_coords.push_back( std::pair<uint,uint>( cur_exon.left(), cur_exon.right()));
   }
   const GenomicFeature & cur_exon = exons[distance(start_pos.cbegin(),up)];
   exon_coords.push_back(  std::pair<uint,uint>( cur_exon.left(), cur_exon.right()));

   auto c = _coords.cbegin();
   ++c;
   uint i = 1;

   while(i< exon_coords.size()-1){
      if(exon_coords[i].first < c->first){
         idx.push_back(i);
         ++i;
      }
      else if(exon_coords[i].first == c->first){
         ++i;
         ++c;
      }
      else{
         assert(false);
      }
   }
   return idx;
}

inline int ExonBin::left_exon_len() const
{
   return _coords.cbegin()->second - _coords.cbegin()->first + 1;
}


inline int ExonBin::effective_len(const std::vector<uint> & seg_lens,
                           const std::vector<uint>& implicit_idx,
                           const int fl,
                           const int rl) const
{

   int gap = fl - 2*rl;

   if(seg_lens.size() == 1){ // single exon
      return seg_lens[0] - fl + 1;
   }
   else if(seg_lens.size() == 2){ //double seg_lens
      return no_gap_ef(seg_lens[0], seg_lens[1], 0, fl);
   }


   else if(seg_lens.size() == 3){ // triple seg_lens
      if(implicit_idx.size() == 1){
         return gap_ef(seg_lens[0], seg_lens[2], seg_lens[1], rl, gap);
      }

      else if(implicit_idx.size() == 0){
         return ( no_gap_ef(seg_lens[0], seg_lens[2], seg_lens[1], fl) - gap_ef(seg_lens[0], seg_lens[2], seg_lens[1], rl, gap));
      }

      else{
         assert(false);
      }
   }
   else if(seg_lens.size() == 4){ // quadruple seg_lens
      int hit14 = gap_ef(seg_lens[0], seg_lens[3], seg_lens[2] + seg_lens[1], rl, gap);
      int hit24 = gap_ef(seg_lens[3], seg_lens[1], seg_lens[2], rl ,gap);
      int hit124 = gap_ef(seg_lens[0]+seg_lens[1], seg_lens[3], seg_lens[2], rl, gap);
      int hit13 = gap_ef(seg_lens[0], seg_lens[2], seg_lens[1], rl, gap);
      int hit134 = gap_ef(seg_lens[0], seg_lens[2]+seg_lens[3], seg_lens[1], rl, gap);

      if(implicit_idx.size() == 0){
         int hit_all_124 = hit124 - hit14 - hit24;
         int hit_all_134 = hit134 - hit14 - hit13;
         int total = no_gap_ef(seg_lens[0], seg_lens[3], seg_lens[1]+seg_lens[2], fl);
         return (total - hit_all_124 - hit_all_134 - hit14);
      }
      else if(implicit_idx.size() == 2){
         return hit14;
      }
      else{
         if(implicit_idx[0] == 1){
            int hit_all_134 = hit134 - hit14 - hit13;
            return hit_all_134;
         }
         else{
            int hit_all_124 = hit124 - hit14 - hit24;
            return hit_all_124;
         }
      }
   }
   else{ // more than quadruple seg_lens
      uint num_inners = seg_lens.size() - 2;
      uint num_pos = 0;
      uint target = pow(2.0, seg_lens.size()) -1;
      for(uint idx : implicit_idx){
         target &= ~(1u << idx);
      }
      for( int i = 1; i != seg_lens[0]+1; ++i){
         uint hit = 1;
         int bp_last = fl - i - accumulate(seg_lens.begin()+1, seg_lens.end()-1,0);
         if(bp_last > *seg_lens.rbegin()) continue;
         if(bp_last < 0) assert(false);
         if(bp_last == 0) break;

         hit |= (1u << (seg_lens.size() -1));//
         //right-end cover
         int last_rest_bp = rl - bp_last;
         uint j = num_inners;
         while(last_rest_bp > 0 && j > 0){
            hit |= (1u << j);
            last_rest_bp = last_rest_bp - seg_lens[j];
            j = j-1;
         }

         //left-end cover
         int first_rest_bp = rl - i;
         j = 1;
         while(first_rest_bp > 0 && j <= num_inners){
            hit |= (1u << j);
            first_rest_bp = first_rest_bp - seg_lens[j];
            j = j+1;
         }

         if(hit == target) num_pos++;
      }
//#ifdef DEBUG
//      cout<<"fl: "<<fl<<"  num pos: "<<num_pos<<endl;
//#endif
      return num_pos;
   }
   return 0;
}
//double ExonBin::read_depth() const
/*
 * This function contains bugs.
 */
//{
//   int bin_len = this->bin_len();
//   double bin_cov = 0;
//   for(auto const & kv:_iso_2_frag_lens){
//      for(auto const &frag: kv.second){
//         bin_cov += frag.first * frag.second;
//      }
//   }
//
//   return bin_cov/bin_len;
//}


inline void ExonBin::add_frag_len(const int iso, const int frag_len, const float mass)
{
   add_read_mass(mass);
   auto ret = _iso_2_frag_lens.find(iso);
   if(ret == _iso_2_frag_lens.end()){
      std::vector<std::pair<int,float>> vec_lens = {std::pair<int, float>(frag_len, mass) };
      _iso_2_frag_lens.emplace(iso, vec_lens);
   }
   else{
      ret->second.push_back(std::pair<int,float>(frag_len, mass));
   }
}

#endif //STRAWBERRY_ISOFORM_H
