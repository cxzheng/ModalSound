/******************************************************************************
 *  File: geo2d.h
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
#ifndef GEOMETRY_GEO2D_H
#   define GEOMETRY_GEO2D_H

namespace geo2d
{

template <typename T>
struct Point2
{
    static const T EPS;

    T   x;
    T   y;

    Point2():x(0), y(0) { }
    Point2(T x, T y):x(x), y(y) { }

    static inline bool is_zero(T v)
    {
        return (v < 0 ? (-v) : (v)) < EPS;
    }
};

template<>
const double Point2<double>::EPS = 1e-16;

template<>
const float Point2<float>::EPS = 1e-8;

/*!
 * compute (P1-P0) x (P2-P0)
 */
template <typename T>
inline T xmult(const Point2<T>& p1, const Point2<T>& p2, const Point2<T>& p0)
{
    return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}

/*!
 * check if the point is on the given segment, including the ending points
 * of the segment
 */
template <typename T>
inline bool point_online_in(const Point2<T>& p, const Point2<T>& l1, const Point2<T>& l2)
{
    return Point2<T>::is_zero(xmult(p,l1,l2)) && 
           (l1.x-p.x)*(l2.x-p.x)<Point2<T>::EPS && 
           (l1.y-p.y)*(l2.y-p.y)<Point2<T>::EPS;
}

template <typename T>
inline bool is_parallel(
        const Point2<T>& u1, const Point2<T>& u2, 
        const Point2<T>& v1, const Point2<T>& v2)
{
    return Point2<T>::is_zero((u1.x-u2.x)*(v1.y-v2.y)-(v1.x-v2.x)*(u1.y-u2.y));
}

template <typename T>
inline bool points_inline(const Point2<T>& p1,
        const Point2<T>& p2, const Point2<T>& p3)
{
    return Point2<T>::is_zero(xmult(p1,p2,p3));
}

template <typename T>
inline bool same_side(
        const Point2<T>& p1, const Point2<T>& p2, 
        const Point2<T>& l1, const Point2<T>& l2)
{
    return xmult(l1,p1,l2) * xmult(l1,p2,l2)>Point2<T>::EPS;
}

/*!
 * Check if two segments intersects with each other. 
 * NOTE: It still returns true if the intersection point is on one of the ending points
 *       of the given segments.
 */
template <typename T>
inline bool seg_intersect_in(
        const Point2<T>& u1, const Point2<T>& u2, 
        const Point2<T>& v1, const Point2<T>& v2)
{
     if (!points_inline(u1,u2,v1) || !points_inline(u1,u2,v2))
         return !same_side(u1,u2,v1,v2) && !same_side(v1,v2,u1,u2);

    return point_online_in(u1,v1,v2) ||
           point_online_in(u2,v1,v2) ||
           point_online_in(v1,u1,u2) ||
           point_online_in(v2,u1,u2);
}

/*!
 * Compute the intersection point between two points.
 */
template <typename T>
inline bool seg_intersect_in(
        const Point2<T>& u1, const Point2<T>& u2, 
        const Point2<T>& v1, const Point2<T>& v2,
        Point2<T>& ret)
{
    if ( is_parallel(u1, u2, v1, v2) ) return false;
    if ( !seg_intersect_in(u1, u2, v1, v2) ) return false;

    ret=u1;
    T t = ((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x)) / 
            ((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
    ret.x += (u2.x-u1.x)*t;
    ret.y += (u2.y-u1.y)*t;
    
    return true;
}

}
#endif
