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

#endif //STRAWBERRY_UTILS_H
