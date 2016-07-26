#ifndef SC_VECTOR_2_HPP
#   define SC_VECTOR_2_HPP

#include <cmath>
#include <iostream>
#include <assert.h>
#include "Tuple2.hpp"

//! Class for three dimensional vector.
template <typename T> 
class Vector2 : public Tuple2<T>
{
    public:
        using Tuple2<T>::x;
        using Tuple2<T>::y;

        /*! Creates and sets to (0,0,0) */
        Vector2() {}
        /*! Creates and sets to (x,y,z) */
        Vector2(T nx, T ny):Tuple2<T>(nx, ny) {}
        /*! Copy constructor. */
        Vector2(const Tuple2<T>& src):Tuple2<T>(src) {}
        Vector2(const Vector2<T>& src):Tuple2<T>(src.x, src.y) {}
        /*! Copy casting constructor. */
        template <typename FromT>
        Vector2(const Tuple2<FromT>& src):Tuple2<T>(src) {}
       
        Vector2<T>& operator = (const Tuple2<T>& rhs) 
        {
            x = rhs.x;
            y = rhs.y;
            return *this;
        }

        Vector2<T> operator * (T rhs) const
        {
            return Vector2<T>(x*rhs, y*rhs);
        }

        Vector2<T> operator / (T rhs) const
        {
            return Vector2<T>(x/rhs, y/rhs);
        }

        Vector2<T> operator + (const Vector2<T>& rhs) const 
        {
            return Vector2<T>(x + rhs.x, y + rhs.y);
        }

        /* this x rhs */
        T cross(const Vector2<T>& rhs) const
        {   return x*rhs.y - rhs.x*y; }

        /*! Get lenght of vector.*/
        T length() const 
        {   return (T)std::sqrt(x*x + y*y); }

        T length_sqr() const 
        {   return x*x + y*y; }

        /*! Normalize vector */
        void normalize() 
        {
            T s = length();
            if ( s > 0 )
            {
                s = (T)1 / s;
                x *= s;
                y *= s;
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
            }
            return ret;
        }

        /*! Dot product of two vectors. */
        T dot(const Vector2<T>& rhs) const 
        {   return x * rhs.x + y * rhs.y; }
};

#endif
