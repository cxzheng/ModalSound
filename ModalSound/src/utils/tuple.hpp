/*
 * =====================================================================================
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * -------------------------------------------------------------------------------------
 *
 *       Filename:  tuple.h
 *
 *    Description:  utilities for tuple
 *
 *        Version:  1.0
 *        Created:  07/05/12 00:02:32
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef UTILS_TUPLE_INC
#   define UTILS_TUPLE_INC

#include <algorithm>
#include "sc/Tuple3.hpp"

template <typename T>
inline void sort_triple(T& a, T& b, T& c)
{
    if ( b < a ) std::swap(a, b);
    if ( c < a ) std::swap(a, c);
    if ( c < b ) std::swap(b, c);
}

template <typename T>
inline int max_in_triple(const T* a)
{
    return a[0] > a[1] ? (a[0] > a[2] ? 0 : 2) : 
                         (a[1] > a[2] ? 1 : 2);
}

template <typename T>
inline int min_in_triple(const T* a)
{
    return a[0] < a[1] ? (a[0] < a[2] ? 0 : 2) :
                         (a[1] < a[2] ? 1 : 2);
}

template <typename T>
inline int idx_in_triple(T a, const Tuple3<T>& c)
{
    return a == c.x ? 0 : 
          (a == c.y ? 1 : 
          (a == c.z ? 2 : -1));
}
#endif

