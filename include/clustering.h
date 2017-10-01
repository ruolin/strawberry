//
// Created by ruolin on 9/27/17.
//

#ifndef STRAWBERRY_CLUSTERING_HPP_H
#define STRAWBERRY_CLUSTERING_HPP_H

#include <vector>
#include <stdexcept>
#include <tuple>

inline size_t one_d_binary_clustering(const std::vector<int>& vec) {
   int total_c1 = 0;
   int total_c2 = 0;
   for (auto const& i : vec) {
      if (i == 0) {
         total_c1++;
      }
      else if (i == 1) {
         total_c2 ++;
      } else {
         throw std::runtime_error("binary clustering must contain only 0 or 1.");
      }
   }
   int l_c1 = 0;
   int l_c2 = 0;
   int r_c1 = total_c1;
   int r_c2 = total_c2;
   int best_score = 0;
   size_t best_idx = 0;
   for (size_t i = 0; i < vec.size(); ++i) {
      if (vec[i] == 0) {
         l_c1 ++;
         r_c1 --;
      } else {
         l_c2 ++;
         r_c2 --;
      }
      int s = std::max(l_c1, l_c2) + std::max(r_c1, r_c2);
      if (s >= best_score) {
         best_idx = i;
         best_score = s;
      }
   }
   return best_idx;
}

#endif //STRAWBERRY_CLUSTERING_HPP_H
