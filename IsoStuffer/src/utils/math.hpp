/******************************************************************************
 *  File: math.hpp
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
/*
 * math.hpp
 * author: Changxi Zheng (cxzheng@cs.cornell.edu)
 */
#ifndef MATH_UTILS_H
#   define MATH_UTILS_H

#ifndef M_PI
#   define M_PI        3.14159265358979323846264338327950288   /* pi */
#endif

#ifndef M_1_PI
#   define M_1_PI      0.318309886183790671537767526745028724  /* 1/pi */
#endif

template <typename T> inline T M_NEG(T a)
{   return -a; }

template <typename T> inline T M_MAX(const T a, const T b)
{   return a > b ? a : b; }

template <typename T> inline T M_MIN(const T a, const T b)
{   return a < b ? a : b; }

template <typename T> inline T M_DEG2RAD(T x)
{   return x * M_PI / 180.0; }

template <typename T> inline T M_RAD2DEG(T x)
{   return x * 180.0 * M_1_PI; }

template <typename T> inline T M_ABS(T a)
{   return a > (T)0 ? a : -a; }

template <typename T> inline T M_TRI(T a)
{   return a * a * a; }

template <typename T> inline T M_SQR(T x)
{   return x*x; }

template <typename T> inline int M_SIGN(T x)
{   return x > (T)0 ? 1 : (x < (T)0 ? -1 : 0); }

inline int gcd(int a, int b)
{
    for(int c;b;c=a,a=b,b=c%b) ;
    return a;
}

inline int lcm(int a, int b)
{
    return a/gcd(a,b)*b;
}

template <typename T>
inline T clamp(T a, T minv, T maxv)
{
    return a <= minv ? minv : (a > maxv ? maxv : a);
}

template <typename T>
inline T gaussian(T v, T mu, T d)
{
    return exp(-M_SQR((v - mu) / d));
}

#endif
