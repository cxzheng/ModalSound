/******************************************************************************
 *  File: Point3.hpp
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
#ifndef GEOMETRY_POINT_3_HPP
#   define GEOMETRY_POINT_3_HPP

#include <assert.h>
#include <cmath>
#include "linearalgebra/Tuple3.hpp"
#include "linearalgebra/Vector3.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

//! A 3-dimensional Point
template <typename T>
class Point3 : public Tuple3<T>
{
    public:
        using Tuple3<T>::x;
        using Tuple3<T>::y;
        using Tuple3<T>::z;

        static const Point3<T>  ORIGIN;

        // ========== Constructors ===========
        Point3() {}
        Point3(T nx, T ny, T nz):Tuple3<T>(nx, ny, nz) {}
        Point3(const Tuple3<T>& src):Tuple3<T>(src) {}

        template <typename FromT>
        Point3(const Tuple3<FromT>& src):Tuple3<T>(src) {}

        Point3<T>& operator = (const Tuple3<T>& rhs) 
        {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            return *this;
        }

        /*! Copy casting operator. */
        template <typename FromT>
        Point3<T>& operator = (const Tuple3<FromT>& rhs)
        {
            x = static_cast<T>(rhs.x);
            y = static_cast<T>(rhs.y);
            z = static_cast<T>(rhs.z);
            return *this;
        }

        /*! distance between two vectors */
        T distance(const Point3<T>& v) const
        {
            T dx = x - v.x;
            T dy = y - v.y;
            T dz = z - v.z;
            return (T)std::sqrt(dx*dx + dy*dy + dz*dz);
        }
        
        T distanceSqr(const Point3<T>& v) const
        {
            T dx = x - v.x;
            T dy = y - v.y;
            T dz = z - v.z;
            return dx*dx + dy*dy + dz*dz;
        }

        /*! Substraction operator */
        Vector3<T> operator - (const Point3<T>& rhs) const 
        {
            return Vector3<T>(x - rhs.x, y - rhs.y, z - rhs.z);
        }

        /*! Addition operator */
        Point3<T> operator + (const Vector3<T>& rhs) const 
        {
            return Point3<T> (x + rhs.x, y + rhs.y, z + rhs.z);
        }

        Point3<T>& operator += (const Vector3<T>& rhs) 
        {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        Point3<T>& operator += (T rhs) 
        {
            x += rhs;
            y += rhs;
            z += rhs;
            return *this;
        }

        Point3<T>& operator *= (T rhs) 
        {
            x *= rhs;
            y *= rhs;
            z *= rhs;
            return *this;
        }

        Point3<T> operator - (T rhs) const 
        {
            return Point3<T>(x - rhs, y - rhs, z - rhs);
        }

        Point3<T> operator + (T rhs) const 
        {
            return Point3<T>(x + rhs, y + rhs, z + rhs);
        }

        Point3<T> operator - () const
        {
            return Point3<T>(-x, -y, -z);
        }

        Point3<T> operator * (T rhs) const
        {
            return Point3<T>(x*rhs, y*rhs, z*rhs);
        }
};

typedef class Point3<float>     Point3f;
typedef class Point3<double>    Point3d;
typedef class Point3<int>       Point3i;

template <typename T> 
const Point3<T> Point3<T>::ORIGIN(0, 0, 0);

#ifdef USE_NAMESPACE
}
#endif

#endif
