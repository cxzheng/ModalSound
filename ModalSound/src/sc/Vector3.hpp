/******************************************************************************
 *  File: Vector3.hpp
 *  A vector in 3D space
 *  Copyright (c) 2007 by Changxi Zheng
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
#ifndef LINEARALGEBRA_VECTOR_3_HPP
#   define LINEARALGEBRA_VECTOR_3_HPP

#include <cmath>
#include <iostream>
#include <assert.h>
#include "Tuple3.hpp"
#include "utils/math.hpp"

#ifdef USE_NAMESPACE
namespace sploosh
{
#endif
    
//! Class for three dimensional vector.
template <typename T> 
class Vector3 : public Tuple3<T>
{
    public:
        using Tuple3<T>::x;
        using Tuple3<T>::y;
        using Tuple3<T>::z;

        static const Vector3<T>  ZERO;

        /*! Creates and sets to (0,0,0) */
        Vector3() {}
        /*! Creates and sets to (x,y,z) */
        Vector3(T nx, T ny, T nz):Tuple3<T>(nx, ny, nz) {}
        /*! Copy constructor. */
        Vector3(const Tuple3<T>& src):Tuple3<T>(src) {}
        /*! Copy casting constructor. */
        template <typename FromT>
        Vector3(const Tuple3<FromT>& src):Tuple3<T>(src) {}
       
        Vector3<T>& operator = (const Tuple3<T>& rhs) 
        {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            return *this;
        }

        /*! Copy casting operator. */
        /*
        template <typename FromT>
        Vector3<T>& operator = (const Tuple3<FromT>& rhs)
        {
            x = static_cast<T>(rhs.x);
            y = static_cast<T>(rhs.y);
            z = static_cast<T>(rhs.z);
            return *this;
        }
        */

        Vector3<T>& operator = (const T rhs)
        {
            x = y = z = rhs;
            return *this;
        }
        /*
        template <typename FromT>
        Vector3<T>& operator = (const FromT rhs)
        {
            x = y = z = static_cast<T>(rhs);
            return *this;
        }
        */

        Vector3<T> operator - (const Vector3<T>& rhs) const 
        {
            return Vector3<T>(x - rhs.x, y - rhs.y, z - rhs.z);
        }

        Vector3<T> operator - () const
        {
            return Vector3<T>(-x, -y, -z);
        }

        /*! Multiplication operator */
        Vector3<T> operator * (const T& rhs) const
        {
            return Vector3<T>(x*rhs, y*rhs, z*rhs);
        }

        Vector3<T> operator * (const Vector3<T>& rhs) const 
        {
            return Vector3<T>(x * rhs.x, y * rhs.y, z * rhs.z);
        }

        Vector3<T> operator / (const T& rhs) const
        {
            return Vector3<T>(x/rhs, y/rhs, z/rhs);
        }

        /*! Dot product of two vectors. */
        T dot(const Vector3<T>& rhs) const 
        {
            return x * rhs.x + y * rhs.y + z * rhs.z;
        }

        T dot(const Tuple3<T>& rhs) const 
        {
            return x * rhs.x + y * rhs.y + z * rhs.z;
        }
       
        /*! Cross product opertor */    
        Vector3<T> cross(const Vector3<T>& rhs) const // $$TESTED
        {
            return Vector3<T>(y * rhs.z - rhs.y * z, 
                              z * rhs.x - rhs.z * x, 
                              x * rhs.y - rhs.x * y);
        }
       
        /*! Get lenght of vector.*/
        T length() const 
        {
            return (T)std::sqrt(x * x + y * y + z * z);
        }
       
        /*!
         * Return square of length.
         * @return length ^ 2
         * @note This method is faster then length(). For comparison
         * of length of two vector can be used just this value, instead
         * of computionaly more expensive length() method.
         */
        T length_sqr() const 
        {
            return x * x + y * y + z * z;
        }
       
        /*! Normalize vector */
        void normalize() 
        {
            T s = length();
            if ( s > 0 )
            {
                s = (T)1 / s;
                x *= s;
                y *= s;
                z *= s;
            }
        }

        /*!
         * Normalize vector and return its original length
         */
        T normalize2()
        {
            T ret = length();
            if ( ret > 0 )
            {
                T s = (T)1 / ret;
                x *= s;
                y *= s;
                z *= s;
            }
            return ret;
        }
       
        /*!
         * Rotate vector around three axis.
         * @param ax Angle (in degrees) to be rotated around X-axis.
         * @param ay Angle (in degrees) to be rotated around Y-axis.
         * @param az Angle (in degrees) to be rotated around Z-axis.
         */
        void rotate(T ax, T ay, T az) 
        {
            T a = (T)cos(DEG2RAD(ax));
            T b = (T)sin(DEG2RAD(ax));
            T c = (T)cos(DEG2RAD(ay));
            T d = (T)sin(DEG2RAD(ay));
            T e = (T)cos(DEG2RAD(az));
            T f = (T)sin(DEG2RAD(az));
            T nx = c * e * x - c * f * y + d * z;
            T ny = (a * f + b * d * e) * x + (a * e - b * d * f) * y - b * c * z;
            T nz = (b * f - a * d * e) * x + (a * d * f + b * e) * y + a * c * z;
            x = nx; y = ny; z = nz;
        }
       
}; // end of Vector3

typedef class Vector3<float>    Vector3f;
typedef class Vector3<double>   Vector3d;
typedef class Vector3<int>      Vector3i;

template <typename T>
const Vector3<T> Vector3<T>::ZERO(0, 0, 0);

#ifdef USE_NAMESPACE
}
#endif

#endif
