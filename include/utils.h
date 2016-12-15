//
// Created by ruolin on 12/10/16.
//

#ifndef STRAWBERRY_UTILS_H
#define STRAWBERRY_UTILS_H

#include<vector>
#include<iostream>

template < class T >
inline std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
{
   typedef typename std::vector<T>::const_iterator const_iterator;
   os << "[";
   for (const_iterator ii = v.begin(); ii != v.end(); ++ii)
   {
      os << " " << *ii;
   }
   os << " ]";
   return os;
}

template<typename Item>
int UniqPushAndReturnIdx(Item& item, std::vector<Item>& container) {
   int idx;
   auto search = std::find(container.begin(), container.end(), item);
   if (search == container.end()) {
      idx = container.size();
      item.id() = idx;
      container.push_back(item);
   } else {
      idx = std::distance(container.begin(), search);
   }
   return idx;
}

template<typename Item>
int PushAndReturnIdx(Item& item, std::vector<Item>& container) {
   int idx = container.size();
   item.id() = idx;
   container.push_back(item);
   return idx;
}

#endif //STRAWBERRY_UTILS_H
