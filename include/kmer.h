//
// Created by Ruolin Liu on 4/15/17.
//

#ifndef STRAWBERRY_KMER_H
#define STRAWBERRY_KMER_H

#include <vector>
#include <iterator>
#include <cmath>
#include <algorithm>

template<typename Seq>
class Kmer {
public:
   static decltype(auto) SortedKmer(const Seq& seq, const int k) {
      assert( k < 32 && k > 0);
      std::vector<uint64_t> sorted_kmers;
      sorted_kmers.reserve(seq.size());
      assert(seq.size() > k);
      auto itr = std::begin(seq);

      //first kmer
      for (int i = 0; i < k; ++i) {
         if (i == 0) {
            sorted_kmers.push_back(ToDna(*itr++));
         } else {
            sorted_kmers.back() <<= 2;
            sorted_kmers.back() |= ToDna(*itr++);
         }
      }

      uint64_t mask = std::pow(2, k*2) - 1;
      //rest
      for (; itr != std::end(seq); ++itr) {
         uint64_t tmp = sorted_kmers.back();
         tmp <<= 2;
         tmp |= ToDna(*itr);
         tmp &= mask;
         sorted_kmers.push_back(tmp);
      }
      std::sort(sorted_kmers.begin(), sorted_kmers.end());
      return sorted_kmers;
   }

   static double Entropy(Seq const& seq, int k) {
      /*Kmer entropy*/
      auto sorted_kmer = SortedKmer(seq, k);
      size_t total = sorted_kmer.size();
      double counter = 1.0;
      double sum = 0.0;
      for (size_t i = 1; i < total; ++i) {
         if (sorted_kmer[i] != sorted_kmer[i-1]) {
            double p = counter/total;
            sum -= p * log(p);
            counter = 1.0;
         } else {
            counter += 1.0;
         }
      }
      double p = counter/total;
      sum -= p * log(p);
      return sum;
   }

   template<typename Itr>
   static double GCRatio(Itr b, Itr e){
      assert (b != e);
      int gc_counter = 0;
      int total = 0;
      for (; b != e; ++b) {
         ++total;
         gc_counter += ToDna2(*b);
      }
      return (double) gc_counter / total;
   }

   template<typename Itr>
   static bool HighGCStrech(Itr b, Itr e, int w , double cutoff) {
      assert(cutoff <= 1.0);
      auto len = std::distance(b, e);
      assert ( w < len);
      while (b + w <= e) {
         if (GCRatio(b, b + w) > cutoff) return true;
         ++b;
      }
      return false;
   }

private:
   static decltype(auto) ToDna2(char c) {
      switch (c) {
         case 'C':
         case 'c':
         case 'G':
         case 'g':
         case 1u:
         case 2u:
            return 1u;
         default:
            return 0u;
      }

   }
   static decltype(auto) ToDna(char c) {
      switch (c) {
         case 'A':
         case 'a':
            return 0u;
         case 'C':
         case 'c':
            return 1u;
         case 'G':
         case 'g':
            return 2u;
         case 'T':
         case 't':
            return 3u;
         default:
            return 0u;
      }
   }
};

#endif //STRAWBERRY_KMER_H
