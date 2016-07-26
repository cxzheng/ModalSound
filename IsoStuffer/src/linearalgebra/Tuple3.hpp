/******************************************************************************
 *  File: Tuple3.hpp
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
#ifndef GEOMETRY_TUPLE_3_HPP
#   define GEOMETRY_TUPLE_3_HPP

#include <iostream>
#include <assert.h>

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T>
class Tuple3
{
    public:
        typedef T   element;

        union
        {
            T x;    // First element of vector, alias for X-coordinate.
            T u;    // for velocity
            T r;
        };

        union
        {
            T y; 
            T v;
            T theta;    // for spheric coordinates
            T g;
        };

        union
        {
            T z;
            T w;
            T phi;
            T b;
        };

        // ========== Constructors ===========
        Tuple3():x(0), y(0), z(0) {}
        Tuple3(T nx, T ny, T nz):x(nx), y(ny), z(nz) {}
        Tuple3(const Tuple3<T>& src):x(src.x), y(src.y), z(src.z) {}

        template <typename FromT>
        Tuple3(const Tuple3<FromT>& src):x(static_cast<T>(src.x)),
                                         y(static_cast<T>(src.y)),
                                         z(static_cast<T>(src.z))
        {}

        /*! Directly set the fields */
        template <typename FromT>
        void set(FromT x, FromT y, FromT z)
        {
            this->x = static_cast<T>(x);
            this->y = static_cast<T>(y);
            this->z = static_cast<T>(z);
        }

        void zero()
        {
            x = static_cast<T>(0);
            y = static_cast<T>(0);
            z = static_cast<T>(0);
        }

        //================ operators ==============
        /*! Copy operator */
        Tuple3<T>& operator = (const Tuple3<T>& rhs) 
        {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            return *this;
        }
       
        /*! Copy casting operator. */
        template <typename FromT>
        Tuple3<T>& operator = (const Tuple3<FromT>& rhs)
        {
            x = static_cast<T>(rhs.x);
            y = static_cast<T>(rhs.y);
            z = static_cast<T>(rhs.z);
            return *this;
        }

        /*!
         * Array access operator
         * @param n Array index
         * @return For n = 0, reference to x coordinate, n = 1
         * reference to y, else reference to z 
         * y coordinate.
         */
        T& operator [] (int n)
        {
            assert(n >= 0 && n <= 2);
            //return n == 0 ? x : (n == 1 ? y : z);
            return ((T*)this)[n];
        }
       
        /*! Addition operator */
        Tuple3<T> operator + (const Tuple3<T>& rhs) const 
        {
            return Tuple3<T> (x + rhs.x, y + rhs.y, z + rhs.z);
        }
       
        /*! Substraction operator */
        Tuple3<T> operator - (const Tuple3<T>& rhs) const 
        {
            return Tuple3<T>(x - rhs.x, y - rhs.y, z - rhs.z);
        }
       
        Tuple3<T> operator - () const
        {
            return Tuple3<T>(-x, -y, -z);
        }

        /*! Multiplication operator */
        Tuple3<T> operator * (const Tuple3<T>& rhs) const 
        {
            return Tuple3<T> (x * rhs.x, y * rhs.y, z * rhs.z);
        }
       
        /*! Division operator */
        Tuple3<T> operator / (const Tuple3<T>& rhs) const 
        {
            return Tuple3<T> (x / rhs.x, y / rhs.y, z / rhs.z);
        }
       
        /*! Addition operator */
        Tuple3<T>& operator += (const Tuple3<T>& rhs) 
        {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }
       
        /*! Substraction operator */
        Tuple3<T>& operator -= (const Tuple3<T>& rhs) 
        {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }
       
        template <typename FromT>
        Tuple3<T>& operator -= (const Tuple3<FromT>& rhs)
        {
            x -= static_cast<T>(rhs.x);
            y -= static_cast<T>(rhs.y);
            z -= static_cast<T>(rhs.z);
            return  *this;
        }

        /*! Multiplication operator */
        Tuple3<T>& operator *= (const Tuple3<T>& rhs) 
        {
            x *= rhs.x;
            y *= rhs.y;
            z *= rhs.z;
            return *this;
        }
       
        /*! Division operator */
        Tuple3<T>& operator /= (const Tuple3<T>& rhs) 
        {
            x /= rhs.x;
            y /= rhs.y;
            z /= rhs.z;
            return *this;
        }

        //============== scalar vector operator ==============
        /*! Addition operator */
        Tuple3<T> operator + (T rhs) const 
        {
            return Tuple3<T> (x + rhs, y + rhs, z + rhs);
        }
       
        /*! Substraction operator */
        Tuple3<T> operator - (T rhs) const 
        {
            return Tuple3<T> (x - rhs, y - rhs, z - rhs);
        }
       
        /*! Multiplication operator */
        Tuple3<T> operator * (T rhs) const 
        {
            return Tuple3<T> (x * rhs, y * rhs, z * rhs);
        }
       
        friend Tuple3<T> operator * (T lhs, const Tuple3<T>& rhs)
        {
            return Tuple3<T>(lhs*rhs.x, lhs*rhs.y, lhs*rhs.z);
        }

        /*! Division operator */
        Tuple3<T> operator / (T rhs) const 
        {
            return Tuple3<T> (x / rhs, y / rhs, z / rhs);
        }
       
        /*! Addition operator */
        Tuple3<T>& operator += (T rhs) 
        {
            x += rhs;
            y += rhs;
            z += rhs;
            return *this;
        }
       
        /*! Substraction operator */
        Tuple3<T>& operator -= (T rhs) 
        {
            x -= rhs;
            y -= rhs;
            z -= rhs;
            return *this;
        }
       
        /*! Multiplication operator */
        Tuple3<T>& operator *= (T rhs) 
        {
            x *= rhs;
            y *= rhs;
            z *= rhs;
            return *this;
        }
       
        /*! Division operator */
        Tuple3<T>& operator /= (T rhs) 
        {
            x /= rhs;
            y /= rhs;
            z /= rhs;
            return *this;
        }

        bool equals(const Tuple3<T>& rhs) const
        {
            return x == rhs.x && y == rhs.y && z == rhs.z;
        }

        /*!
         * this = s*this + t1
         */
        void selfScaleAdd(T s, const Tuple3<T>& t1)
        {
            x = x * s + t1.x;
            y = y * s + t1.y;
            z = z * s + t1.z;
        }

        /*!
         * provide the similar method with j3d.vecmath.Tuple2d
         * this = s*t1 + t2
         */
        void scaleAdd(T s, const Tuple3<T>& t1, const Tuple3<T>& t2)
        {
            x = s*t1.x + t2.x;
            y = s*t1.y + t2.y;
            z = s*t1.z + t2.z;
        }

        /*! this = this + s*t1 */
        void scaleAdd(T s, const Tuple3<T>& t1)
        {
            x += s*t1.x;
            y += s*t1.y;
            z += s*t1.z;
        }

        T norm() const
        {
            return (T)sqrt(x*x + y*y + z*z);
        }

        T normSqr() const 
        {
            return x * x + y * y + z * z;
        }

        void clamp(const Tuple3<T>& minV, const Tuple3<T>& maxV)
        {
            x = x <= minV.x ? minV.x : (x > maxV.x ? maxV.x : x);
            y = y <= minV.y ? minV.y : (y > maxV.y ? maxV.y : y);
            z = z <= minV.z ? minV.z : (z > maxV.z ? maxV.z : z);
        }

        /*!
         * Linear interpolation of two vectors
         * @param fact Factor of interpolation. For translation from position
         *        of this vector to vector r, values of factor goes from 0.0 to 1.0.
         * @param r Second vector for interpolation 
         * @note However values of fact parameter are reasonable only in interval
         *       [0.0 , 1.0], you can pass also values outside of this interval and you
         *       can get result (extrapolation?)
         */
        Tuple3<T> lerp(T fact, const Tuple3<T>& r) const
        {
            return (*this) + (r - (*this)) * fact;    
        }

        /*!
         * Conversion to pointer operator
         * @return Pointer to internaly stored (in managment of class Tuple3<T>)
         *         used for passing Tuple3<T> values to gl*3[fd] functions.
         */
        operator T*() { return (T*)this; }
       
        /*!
         * Conversion to pointer operator
         * @return Constant Pointer to internaly stored (in managment of class Tuple3<T>)
         * used for passing Tuple3<T> values to gl*3[fd] functions.
         */    
        operator const T*() const { return (const T*) this; }
       
        //============== stream operator =================
        /*! Output to stream operator */
        friend std::ostream& operator<<(std::ostream& lhs, const Tuple3<T>& rhs) 
        {
            lhs << "[" << rhs.x << "," << rhs.y << "," << rhs.z  << "]";
            return lhs;
        }
       
        /*! Input from stream */
        friend std::istream& operator>>(std::istream& in, Tuple3<T>& obj)
        {
            in >> obj.x >> obj.y >> obj.z;
            return in;
        }
};

typedef class Tuple3<int>           Tuple3i;
typedef class Tuple3<float>         Tuple3f;
typedef class Tuple3<double>        Tuple3d;
typedef class Tuple3<unsigned int>  Tuple3ui;

#ifdef USE_NAMESPACE
}
#endif

#endif
