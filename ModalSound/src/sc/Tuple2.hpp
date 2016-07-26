#ifndef SC_TUPLE_2_INC
#   define SC_TUPLE_2_INC

#include <iostream>
#include <assert.h>

#ifdef USE_NAMESPACE
namespace sploosh
{
#endif

template <typename T>
class Tuple2
{
    public:
        typedef T   element;

        union
        {
            T x;    // First element of vector, alias for X-coordinate.
            T u;    // for velocity
            T c;    // column
        };

        union
        {
            T y; 
            T v;
            T r;    // row
        };

        // ========== Constructors ===========
        Tuple2():x(0), y(0) {}
        Tuple2(T nx, T ny):x(nx), y(ny) {}
        Tuple2(const Tuple2<T>& src):x(src.x), y(src.y) {}

        template <typename FromT>
        Tuple2(const Tuple2<FromT>& src):x(static_cast<T>(src.x)),
                                         y(static_cast<T>(src.y))
        {}

        /*! Directly set the fields */
        template <typename FromT>
        void set(FromT x, FromT y)
        {
            this->x = static_cast<T>(x);
            this->y = static_cast<T>(y);
        }

        void zero()
        {
            x = static_cast<T>(0);
            y = static_cast<T>(0);
        }

        //================ operators ==============
        /*! Copy operator */
        Tuple2<T>& operator = (const Tuple2<T>& rhs) 
        {
            x = rhs.x;
            y = rhs.y;
            return *this;
        }

        Tuple2<T>& operator *= (T rhs) 
        {
            x *= rhs;
            y *= rhs;
            return *this;
        }

        /*! Addition operator */
        Tuple2<T>& operator += (T rhs) 
        {
            x += rhs;
            y += rhs;
            return *this;
        }
        Tuple2<T>& operator += (const Tuple2<T> rhs) 
        {
            x += rhs.x;
            y += rhs.y;
            return *this;
        }

        /*! Substraction operator */
        Tuple2<T> operator - (const Tuple2<T>& rhs) const 
        {
            return Tuple2<T>(x - rhs.x, y - rhs.y);
        }

        Tuple2<T> operator - () const
        {
            return Tuple2<T>(-x, -y, -z);
        }

        /*! Multiplication operator */
        Tuple2<T> operator * (T rhs) const 
        {   return Tuple2<T>(x * rhs, y * rhs); }

        /*! Division operator */
        Tuple2<T>& operator /= (T rhs) 
        {
            x /= rhs;
            y /= rhs;
            return *this;
        }

        T norm() const
        {
            return (T)sqrt(x*x + y*y);
        }

        /*!
         * Array access operator
         * @param n Array index
         * @return For n = 0, reference to x coordinate, n = 1
         * reference to y coordinate.
         */
        T& operator [] (int n)
        {
            assert(n >= 0 && n <= 1);
            return ((T*)this)[n];
        }
        //============== stream operator =================
        /*! Output to stream operator */
        friend std::ostream& operator<<(std::ostream& lhs, const Tuple2<T>& rhs) 
        {
            lhs << "[" << rhs.x << "," << rhs.y << "]";
            return lhs;
        }
};

typedef class Tuple2<int>           Tuple2i;
typedef class Tuple2<float>         Tuple2f;
typedef class Tuple2<double>        Tuple2d;
typedef class Tuple2<unsigned int>  Tuple2ui;

#ifdef USE_NAMESPACE
}
#endif
#endif

