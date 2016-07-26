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
