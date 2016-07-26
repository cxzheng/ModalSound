/******************************************************************************
 *  File: tuple.hpp
 *
 *  This file is part of isostuffer
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
 ******************************************************************************/
#ifndef UTILS_TUPLE_HPP
#   define UTILS_TUPLE_HPP

#include <algorithm>

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

#endif
