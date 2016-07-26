#ifndef GEOMETRY_POINT_2_HPP
#   define GEOMETRY_POINT_2_HPP

#include <assert.h>
#include <cmath>
#include <complex>
#include "sc/Tuple2.hpp"
#include "sc/Vector2.hpp"

#ifdef USE_NAMESPACE
namespace sploosh
{
#endif

//! A 3-dimensional Point
template <typename T>
class Point2 : public Tuple2<T>
{
    public:
        using Tuple2<T>::x;
        using Tuple2<T>::y;

        static const Point2<T>  ORIGIN;

        // ========== Constructors ===========
        Point2() {}
        Point2(T nx, T ny):Tuple2<T>(nx, ny) {}
        Point2(const Tuple2<T>& src):Tuple2<T>(src) {}

        Point2<T>& operator = (const Tuple2<T>& rhs) 
        {
            x = rhs.x;
            y = rhs.y;
            return *this;
        }

        /*! distance between two points */
        T distance(const Point2<T>& v) const
        {
            T dx = x - v.x;
            T dy = y - v.y;
            return (T)std::sqrt(dx*dx + dy*dy);
        }

        T distance_sqr(const Point2<T>& v) const
        {
            T dx = x - v.x;
            T dy = y - v.y;
            return dx*dx + dy*dy;
        }

        /*! Substraction operator */
        Vector2<T> operator - (const Point2<T>& rhs) const 
        {
            return Vector2<T>(x - rhs.x, y - rhs.y);
        }
        /*! Addition operator */
        Point2<T> operator + (const Vector2<T>& rhs) const 
        {
            return Point2<T>(x + rhs.x, y + rhs.y);
        }

        Point2<T> operator + (const Point2<T>& rhs) const 
        {
            return Point2<T>(x + rhs.x, y + rhs.y);
        }

        Point2<T>& operator *= (T rhs) 
        {
            x *= rhs;
            y *= rhs;
            return *this;
        }

        Point2<T> operator * (T rhs) const
        {
            return Point2<T>(x*rhs, y*rhs);
        }

        Point2<T> operator / (T rhs) const
        {
            return Point2<T>(x/rhs, y/rhs);
        }

        // std::complex<T> toComplex() const
        // {
        //     return std::complex<T> (x,y);
        // }
};

typedef class Point2<float>     Point2f;
typedef class Point2<double>    Point2d;
typedef class Point2<int>       Point2i;

template <typename T> 
const Point2<T> Point2<T>::ORIGIN(0, 0, 0);

#ifdef USE_NAMESPACE
}
#endif

#endif

