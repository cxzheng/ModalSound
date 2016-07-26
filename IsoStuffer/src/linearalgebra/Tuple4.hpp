/******************************************************************************
 *  File: Tuple4.hpp
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
#ifndef VECMATH_VECTOR_H
#   define VECMATH_VECTOR_H

#include <cmath>
#include <iostream>
#include <assert.h>
#include "utils/math.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif
    
#define _epsilon    1e-10
template <typename T>
class Tuple4
{
    public:
        typedef T   element;

        union 
        { 
            T r;        // First element of vector, alias for R-coordinate.
            T x;        // First element of vector, alias for X-coordinate.
        };

        union 
        { 
            T g;        // Second element of vector, alias for G-coordinate.
            T y;        // Second element of vector, alias for Y-coordinate.
        };
 
        union 
        {
            T b;        // Third element of vector, alias for B-coordinate.
            T z;        // Third element of vector, alias for Z-coordinate.
        };
 
        union 
        {
            T a;        // Fourth element of vector, alias for A-coordinate.
            T w;        // First element of vector, alias for W-coordinate.
        };
 
        // ====== constructors ====== 
        Tuple4() : x(0),y(0),z(0),w(0) { }


        //! Creates and sets to (x,y,z,z)
        Tuple4(T nx, T ny, T nz, T nw):x(nx), y(ny), z(nz), w(nw) { }

        //! Copy constructor.
        Tuple4(const Tuple4<T>& src):
            x(src.x), y(src.y), z(src.z), w(src.w) { }

        template <typename FromT>
        void set(FromT x, FromT y, FromT z, FromT w)
        {
            this->x = static_cast<T>(x);
            this->y = static_cast<T>(y);
            this->z = static_cast<T>(z);
            this->w = static_cast<T>(w);
        }

        //! Copy casting constructor.
        template <typename FromT>
        Tuple4(const Tuple4<FromT>& src):
                x(static_cast<T>(src.x)),
                y(static_cast<T>(src.y)),
                z(static_cast<T>(src.z)),
                w(static_cast<T>(src.w)) { }

        // ====== access operators ======
        //! Copy operator
        Tuple4<T>& operator = (const Tuple4<T>& rhs) 
        {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            w = rhs.w;
            return *this;
        }

        //! Copy casting operator
        template<typename FromT>
        Tuple4<T>& operator = (const Tuple4<FromT>& rhs) 
        {
            x = static_cast<T>(rhs.x);
            y = static_cast<T>(rhs.y);
            z = static_cast<T>(rhs.z);
            w = static_cast<T>(rhs.w);
            return *this;
        }

        /*!
         * Array access operator
         * For n = 0, reference to x coordinate, n = 1
         * reference to y coordinate, n = 2 reference to z,  
         * else reference to w coordinate.
         */
        T & operator [] (int n) 
        {
            assert(n >= 0 && n <= 3);
            //return n==0 ? x : 
            //      (n==1 ? y : 
            //      (n==2 ? z : w));
            return ((T*)this)[n];
        }


        // ====== vector aritmetic operator ======
        //! Addition operator
        Tuple4<T> operator + (const Tuple4<T>& rhs) const 
        {
            return Tuple4<T>(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w);
        }

        //! Substraction operator
        Tuple4<T> operator - (const Tuple4<T>& rhs) const 
        {
            return Tuple4<T>(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
        }

        Tuple4<T> operator - () const
        {
            return Tuple4<T>(-x, -y, -z, -w);
        }

        //! Multiplication operator
        Tuple4<T> operator * (const Tuple4<T> rhs) const 
        {
            return Tuple4<T>(x * rhs.x, y * rhs.y, z * rhs.z, w * rhs.w);
        }

        //! Division operator
        Tuple4<T> operator / (const Tuple4<T>& rhs) const 
        {
            return Tuple4<T> (x / rhs.x, y / rhs.y, z / rhs.z, w / rhs.w);
        }

        //! Addition operator
        Tuple4<T>& operator += (const Tuple4<T>& rhs) 
        {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            w += rhs.w;
            return *this;
        }

        //! Substraction operator
        Tuple4<T>& operator -= (const Tuple4<T>& rhs) 
        {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            w -= rhs.w;
            return *this;
        }

        //! Multiplication operator
        Tuple4<T>& operator *= (const Tuple4<T>& rhs) 
        {
            x *= rhs.x;
            y *= rhs.y;
            z *= rhs.z;
            w *= rhs.w;
            return *this;
        }

        //! Division operator
        Tuple4<T>& operator /= (const Tuple4<T>& rhs) 
        {
            x /= rhs.x;
            y /= rhs.y;
            z /= rhs.z;
            w /= rhs.w;
            return *this;
        }

        // ====== equiality operator ======
        //! Equality test operator
        bool operator == (const Tuple4<T>& rhs) const 
        {
            return std::fabs(x - rhs.x) < _epsilon 
                && std::fabs(y - rhs.y) < _epsilon 
                && std::fabs(z - rhs.z) < _epsilon 
                && std::fabs(w - rhs.w) < _epsilon;
        }

        bool equals(const Tuple4<T>& rhs) const
        {
            return x == rhs.x && y == rhs.y && z == rhs.z && w == rhs.w;
        }

        //! Inequality test operator
        bool operator != (const Tuple4<T>& rhs) const 
        { 
            return !(*this == rhs); 
        }

        // ====== scalar vector operator ======
        //! Addition operator
        Tuple4<T> operator + (T rhs) const 
        {
            return Tuple4<T>(x + rhs, y + rhs, z + rhs, w + rhs);
        }

        //! Substraction operator
        Tuple4<T> operator - (T rhs) const 
        {
            return Tuple4<T>(x - rhs, y - rhs, z - rhs, w - rhs);
        }

        //! Multiplication operator
        Tuple4<T> operator * (T rhs) const 
        {
            return Tuple4<T>(x*rhs, y*rhs, z*rhs, w*rhs);
        }

        friend Tuple4<T> operator * (T lhs, const Tuple4<T>& rhs)
        {
            return Tuple4<T>(lhs*rhs.x, lhs*rhs.y, lhs*rhs.z, lhs*rhs.w);
        }

        //! Division operator
        Tuple4<T> operator / (T rhs) const 
        {
            return Tuple4<T>(x / rhs, y / rhs, z / rhs, w / rhs);
        }

        //! Addition operator
        Tuple4<T>& operator += (T rhs) 
        {
            x += rhs;
            y += rhs;
            z += rhs;
            w += rhs;
            return *this;
        }

        //! Substraction operator
        Tuple4<T>& operator -= (T rhs) 
        {
            x -= rhs;
            y -= rhs;
            z -= rhs;
            w -= rhs;
            return * this;
        }

        //! Multiplication operator
        Tuple4<T>& operator *= (T rhs) 
        {
            x *= rhs;
            y *= rhs;
            z *= rhs;
            w *= rhs;
            return *this;
        }

        //! Division operator
        Tuple4<T>& operator /= (T rhs) 
        {
            x /= rhs;
            y /= rhs;
            z /= rhs;
            w /= rhs;
            return *this;
        }

        // ====== size operations ======
        //! Get lenght of vector.
        T length() const
        {
            return (T)std::sqrt(x * x + y * y + z * z + w * w);
        }

        /**
         * Normalize vector
         */
        void normalize() 
        {
            T s = length();
            x /= s;
            y /= s;
            z /= s;
            w /= s;
        }
 
        /**
         * Normalize vector and return its original length
         */
        T normalize2()
        {
            T s = length();
            x /= s;
            y /= s;
            z /= s;
            w /= s;
            return s;
        }

        //! Return square of length.
        T lengthSqr() const 
        {
            return x * x + y * y + z * z + w * w;
        }
 
        /*!
         * Linear interpolation of two vectors
         * \param fact Factor of interpolation. For translation from positon
         *        of this vector to vector r, values of factor goes from 0.0 to 1.0.
         * \param r Second vector for interpolation 
         * \note Hovewer values of fact parameter are reasonable only in interval
         * [0.0 , 1.0], you can pass also values outside of this interval and you
         * can get result (extrapolation?)
         */
        Tuple4<T> lerp(T fact, const Tuple4<T>& r) const
        {
            return (*this) + (r - (*this)) * fact;    
        }
 
        // ====== conversion ======
        //! Conversion to pointer operator
        operator T*() { return (T*)this; }
 
        /*!
         * Conversion to pointer operator
         * return Constant Pointer to internaly stored (in managment of class Tuple4<T>)
         * used for passing Tuple4<T> values to gl*4[fd] functions.
         */
        operator const T* () const { return (const T*)this; }

        //! Output to stream operator
        friend std::ostream& operator<<(std::ostream& lhs, const Tuple4<T>& rhs) 
        {
            lhs << "[" << rhs.x << "," << rhs.y << "," << rhs.z << "," << rhs.w << "]";
            return lhs;
        }
      
        //! Input from stream operator
        friend std::istream& operator>>(std::istream& in, Tuple4<T>& obj)
        {
             in >> obj.x >> obj.y >> obj.z >> obj.w;
             return in;
        }
};

typedef class Tuple4<int>           Tuple4i;
typedef class Tuple4<unsigned int>  Tuple4ui;
typedef class Tuple4<float>         Tuple4f;
typedef class Tuple4<double>        Tuple4d;

#ifdef USE_NAMESPACE
}
#endif

#endif
