/******************************************************************************
 *  File: Point3KDStruct.hpp
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
#ifndef GEOMETRY_POINT3_KDTREE_STRUCT_HPP
#   define GEOMETRY_POINT3_KDTREE_STRUCT_HPP

#include <assert.h>
#include "Point3.hpp"

/*
 * The structures defined here are used when constructing KD-tree, serving as 
 * the input class parameters
 */

template <typename T>
struct Point3Acc
{
    typedef T   result_type;

    inline T operator()(const Point3<T>& v, const size_t N) const
    {
        assert(N < 3);
        return v[N];
    }
};

template <typename T>
struct Point3DistSqr
{
    inline T operator()(const Point3<T>& a, const Point3<T>& b) const
    {   return a.distanceSqr(b); }
};

#endif
