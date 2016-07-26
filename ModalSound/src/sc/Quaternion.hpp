/******************************************************************************
 *  File: Quaternion.hpp
 *  Quaternion in 3D space
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
#ifndef SC_QUATERNION_INC
#   define  SC_QUATERNION_INC

#include <cmath>
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "utils/math.hpp"
#include "generic/precision_type.hpp"

#ifdef USE_NAMESPACE
namespace sploosh
{
#endif

//#define _epsilon    1e-10

template <typename T>
class Quaternion
{
    public:
        typedef T   element;

        T w;                // Real part of Quaternion.
        Vector3<T> v;       // Complex part of Quaternion.
 
        // ====== constructor ======
        //! Create identity Quaternion as default
        Quaternion():w(1), v(0,0,0) { }
        //! Copy constructor
        Quaternion(const Quaternion<T>& q):w(q.w), v(q.v) { }
        //! Copy casting constructor.
        template <typename FromT>
        Quaternion(const Quaternion<FromT>& q):
                w(static_cast<T>(q.w)), v(q.v) { }
        /*!
         * Creates Quaternion object from real part w and complex part v.
         * \param w Real part of Quaternion.
         * \param v Complex part of Quaternion (xi + yj + zk).
         */
        Quaternion(T w, const Vector3<T>& v): w(w), v(v) { }
        //! Creates Quaternion object from value (w + xi + yj + zk).
        Quaternion(T w, T x, T y, T z): w(w), v(x,y,z) { }
        //! Copy operator
        Quaternion<T>& operator = (const Quaternion<T>& rhs)
        {
            v = rhs.v;
            w = rhs.w;
            return *this;
        }
        //! Copy convert operator
        /*
        template<typename FromT>
        Quaternion<T>& operator = (const Quaternion<FromT>& rhs)
        {
            w = static_cast<T>(rhs.w);
            v = rhs.v;
            return *this;
        }
        */

        //! Addition operator
        Quaternion<T> operator + (const Quaternion<T>& rhs) const
        {
            return Quaternion<T>(w + rhs.w, v + rhs.v);  
        }
        //! Multiplication operator
        Quaternion<T> operator * (const Quaternion<T>& rhs) const   // $TESTED$
        {
            return Quaternion<T>(
                    w * rhs.w   - v.x * rhs.v.x - v.y * rhs.v.y - v.z * rhs.v.z,
                    w * rhs.v.x + v.x * rhs.w   + v.y * rhs.v.z - v.z * rhs.v.y,
	            w * rhs.v.y + v.y * rhs.w   + v.z * rhs.v.x - v.x * rhs.v.z,
	            w * rhs.v.z + v.z * rhs.w   + v.x * rhs.v.y - v.y * rhs.v.x);
        }

        //! Multiplication operator
        Quaternion<T> operator * (T rhs) const 
        { 
            return Quaternion<T>(w*rhs, v*rhs);
        }
        //! Substraction operator
        Quaternion<T> operator - (const Quaternion<T>& rhs) const
        {
            return Quaternion<T>(w - rhs.w, v - rhs.v);
        }
        //! Addition operator
        Quaternion<T>& operator += (const Quaternion<T>& rhs)
        {
            w += rhs.w;
            v += rhs.v;
            return *this;
        }
        //! Substraction operator
        Quaternion<T>& operator -= (const Quaternion<T>& rhs)
        {
            w -= rhs.w;
            v -= rhs.v;
            return *this;
        }
        //! Multiplication operator
        Quaternion<T>& operator *= (const Quaternion<T>& rhs)
        {
            (*this) = (*this) * rhs;
            return *this;
        }
        //! Multiplication operator
        Quaternion<T>& operator *= (T rhs)
        {
            w *= rhs;
            v *= rhs;
            return *this;
        }

        /*!  this += s*t1 */
        void scaleAdd(T s, const Quaternion<T>& t1)
        {
            w += s*t1.w;
            v.scaleAdd(s, t1.v);
        }

        //! Equality test operator
        /*!
         * \param rhs Right hand side argument of binary operator.
         * \note Test of equality is based of threshold EPS value. To be two
         *       values equal, must satisfy this condition | lws - rhs | < EPS,
         *       for all Quaternion coordinates.
         */
        bool operator == (const Quaternion<T>& rhs) const
        {
           return std::fabs(w - rhs.w) < PrecisionType<T>::EPS && v == rhs.v;
        }
 
        bool equals(const Quaternion<T>& rhs) const
        {
            return w == rhs.w && v.equals(rhs.v);
        }
        //! Inequality test operator
        bool operator != (const Quaternion<T>& rhs) const 
        { 
            return ! (*this == rhs); 
        }

        //! Get lenght of Quaternion.
        T length() const
        {
            return (T)std::sqrt(w*w + v.length_sqr());    
        }
 
        /*
         * if this quaternion is normalized, then conjugate represents the 
         * inverse rotation of this quaternion
         */
        Quaternion<T> conjugate() const
        {
            return Quaternion<T>(w, -v);
        }

        /*!
         * Return square of length.
         * \return lenght ^ 2
         * \note This method is faster then length(). For comparsion
         * of length of two Quaternion can be used just this value, instead
         * of computionaly more expensive length() method.
         */
        T length_sqr()
        {
            return w*w + v.length_sqr();    
        }
 
        //! Normalize Quaternion
        void normalize()    // $TESTED$
        {
            T len = length_sqr();
            if ( fabs(len - (T)1) > PrecisionType<T>::EPS )
            {
                len = 1. / sqrt(len);
                w *= len;
                v *= len;
            }
        }
 
        /*!
         * Creates Quaternion as rotation around axis.
         * \param axis Unit vector expressing axis of rotation.
         * \param angleDeg Angle of rotation around axis (in degrees).
         */
        static Quaternion<T> from_axis_rot_D(const Vector3<T>& axis, double angleDeg)
        {
            double angleRad = M_DEG2RAD(angleDeg);
            double t = std::sin(angleRad/2);
            double s = std::cos(angleRad/2);
            return Quaternion<T>((T)s, axis*(T)t);
        }
 
        /*!
         * Creates Quaternion as rotation around axis.
         * \param axis Unit vector expressing axis of rotation.
         * \param angleDeg Angle of rotation around axis (in radians).
         */
        static Quaternion<T> from_axis_rot_R(const Vector3<T>& axis, double angle)
        {
            double t = std::sin(angle/2.);
            double s = std::cos(angle/2.);
            return Quaternion<T>((T)s, axis*(T)t);
        }

        /*!
         * Creates Quaternion for eulers angles.
         * \param x Rotation around x axis (in degrees).
         * \param y Rotation around y axis (in degrees).
         * \param z Rotation around z axis (in degrees).
         */
        static Quaternion<T> from_euler_angles_D(T x, T y, T z)
        {
            return from_axis_rot_D(Vector3<T>((T)1,(T)0,(T)0),x) * 
                   from_axis_rot_D(Vector3<T>((T)0,(T)1,(T)0),y) * 
                   from_axis_rot_D(Vector3<T>((T)0,(T)0,(T)1),z);
        }
 
        static Quaternion<T> from_matrix3(const Matrix3<T>& mat)
        {
            const T tr = mat.cols[0].x + mat.cols[1].y + mat.cols[2].z;
            if ( tr > 0 )
            {
                T s = sqrt(tr + (T)1) * 2;
                T invs = (T)1 / s;
                return Quaternion<T>((T)0.25*s, 
                        (mat.cols[1].z - mat.cols[2].y) * invs,
                        (mat.cols[2].x - mat.cols[0].z) * invs,
                        (mat.cols[0].y - mat.cols[1].x) * invs);
            }
            else if ( (mat.cols[0].x > mat.cols[1].y) && (mat.cols[0].x > mat.cols[2].z) )
            {
                T s = sqrt((T)1 + mat.cols[0].x - mat.cols[1].y - mat.cols[2].z) *  2;
                T invs = (T)1 / s;
                return Quaternion<T>(
                        (mat.cols[1].z - mat.cols[2].y) * invs,
                        (T)0.25 * s,
                        (mat.cols[0].y + mat.cols[1].x) * invs,
                        (mat.cols[2].x + mat.cols[0].z) * invs);
            }
            else if ( mat.cols[1].y > mat.cols[2].z )
            {
                T s = sqrt((T)1 + mat.cols[1].y - mat.cols[0].x - mat.cols[2].z) * 2;
                T invs = (T)1 / s;
                return Quaternion<T>(
                        (mat.cols[2].x - mat.cols[0].z) * invs,
                        (mat.cols[0].y + mat.cols[1].x) * invs,
                        (T)0.25 * s,
                        (mat.cols[1].z + mat.cols[2].y) * invs);
            }
            else
            {
                T s = sqrt((T)1 + mat.cols[2].z - mat.cols[0].x - mat.cols[1].y) * 2;
                T invs = (T)1 / s;
                return Quaternion<T>(
                        (mat.cols[0].y - mat.cols[1].x) * invs,
                        (mat.cols[2].x + mat.cols[0].z) * invs,
                        (mat.cols[1].z + mat.cols[2].y) * invs,
                        (T)0.25 * s);

            }
        }

        /*!
         * Converts quaternion into rotation matrix.
         * \return Rotation matrix expressing this quaternion.
         * \note the quaternion must bee normalized
         */
        Matrix3<T> to_matrix3() const
        {
            /*
             * ret.at(0,0) = 1 - 2*v.y*v.y - 2*v.z*v.z;
             * ret.at(1,0) = 2*v.x*v.y - 2*w*v.z;
             * ret.at(2,0) = 2*v.x*v.z - 2*w*v.y;
             *
             * ret.at(0,1) = 2*v.x*v.y + 2*w*v.z;
             * ret.at(1,1) = 1 - 2*v.x*v.x - 2*v.z*v.z;
             * ret.at(2,1) = 2*v.y*v.z - 2*w*v.x;
             *
             * ret.at(0,2) = 2*v.x*v.z - 2*w*v.y;
             * ret.at(1,2) = 2*v.y*v.z + 2*w*v.x;
             * ret.at(2,2) = 1 - 2*v.x*v.x - 2*v.y*v.y;
             */
            T xx = v.x * v.x;
            T xy = v.x * v.y;
            T xz = v.x * v.z;
            T xw = v.x * w;
            
            T yy = v.y * v.y;
            T yz = v.y * v.z;
            T yw = v.y * w;
            
            T zz = v.z * v.z;
            T zw = v.z * w;
            
            return Matrix3<T>(
                    1 - 2*yy - 2*zz, 2*xy - 2*zw, 2*xz + 2*yw,
                    2*xy + 2*zw, 1 - 2*xx - 2*zz, 2*yz - 2*xw,
                    2*xz - 2*yw, 2*yz + 2*xw, 1 - 2*xx - 2*yy);
        }
        /*!
         * Creates Quaternion for eulers angles.
         * \param x Rotation around x axis (in radians).
         * \param y Rotation around y axis (in radians).
         * \param z Rotation around z axis (in radians).
         */
        static Quaternion<T> from_euler_angles_R(T x, T y, T z)
        {
            return fromAxisRotR(Vector3<T>((T)1,(T)0,(T)0),x) * 
                   fromAxisRotR(Vector3<T>((T)0,(T)1,(T)0),y) * 
                   fromAxisRotR(Vector3<T>((T)0,(T)0,(T)1),z);
        }

        T to_axis_rot_D(Vector3<T>& rot) const
        {
            T rad = (T)(2.0 * std::acos(w));
            if ( rad < PrecisionType<T>::EPS )
                rot.set(0, 1, 0);
            else
            {
                rot = v;
                rot.normalize();
            }
            return M_RAD2DEG(rad);
        }

        T to_axis_rot_R(Vector3<T>& rot) const
        {
            T rad = (T)(2.0 * std::acos(w));
            if ( rad < PrecisionType<T>::EPS )
                rot.set(0, 1, 0);
            else
            {
                rot = v;
                rot.normalize();
            }
            return rad;
        }
        //! rotate the given vector
        /*  ret = (s^2-v^2)vec + 2s(v x vec) + 2(v . vec)v
         */
        void rotate_vector(Vector3<T>& vec) const   // $TESTED
        {
            vec = (w*w - v.length_sqr())*vec + 
                  static_cast<T>(2)*w*v.cross(vec) +
                  static_cast<T>(2)*v.dot(vec)*v;
        }

        Vector3<T> rotate(const Vector3<T>& vec) const
        {
            return (w*w - v.length_sqr())*vec + 
                  static_cast<T>(2)*w*v.cross(vec) +
                  static_cast<T>(2)*v.dot(vec)*v;
        }

        /*!
         * Linear interpolation of two Quaternions
         * \param fact Factor of interpolation. For translation from positon
         *        of this vector to Quaternion rhs, values of factor goes from 0.0 to 1.0.
         * \param rhs Second Quaternion for interpolation 
         * \note Values of fact parameter are reasonable only in interval
         * [0.0 , 1.0], you can pass also values outside of this interval and you
         * can get result (extrapolation?)
         */
        Quaternion<T> lerp(T fact, const Quaternion<T>& rhs) const
        {
            return Quaternion<T>((1-fact)*w + fact*rhs.w, v.lerp(fact, rhs.v));
        }

        //! output to standard output stream.
        friend std::ostream& operator << (std::ostream& oss, const Quaternion<T>& q)
        {
            oss << "Re: " << q.w << " Im: " << q.v;
            return oss;
        }

        /*!
         * Computes spherical interpolation between Quaternions (this, q2)
         * using coefficient of interpolation r (in [0, 1]).
         *
         * \param r The ratio of interpolation form this (r = 0) to q2 (r = 1).
         * \param q2 Second Quaternion for interpolation.
         * \return Result of interpolation.
         */
        Quaternion<T> slerp(T r, const Quaternion<T>& q2) const 
        {
            Quaternion<T> ret;
            T cosTheta = w * q2.w + v.x * q2.v.x + v.y *q2.v.y + v.z * q2.v.z;
            if ( fabs(cosTheta) >= 1.0 ) return *this;

            T theta = (T) acos(cosTheta);
            T sinTheta = (T)sqrt(1.0 - cosTheta * cosTheta);
            if (fabs(sinTheta) < PrecisionType<T>::EPS)
            {
                ret.w = (T)(0.5 * w + 0.5 * q2.w);
                ret.v = v.lerp(0.5, q2.v);
            }
            else
            {
                T rA = (T)sin((1.0 - r) * theta) / sinTheta;
                T rB = (T)sin(r * theta) / sinTheta;
                                          
                ret.w = w * rA + q2.w * rB;
                ret.v.x = v.x * rA + q2.v.x * rB;
                ret.v.y = v.y * rA + q2.v.y * rB;
                ret.v.z = v.z * rA + q2.v.z * rB;
            }
            return ret;
        }
};

typedef class Quaternion<float>     Quat4f;
typedef class Quaternion<double>    Quat4d;

#ifdef USE_NAMESPACE
}
#endif

#endif
