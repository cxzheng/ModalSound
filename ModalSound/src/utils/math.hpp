/*
 * math.hpp
 * author: Changxi Zheng (cxzheng@cs.cornell.edu)
 */
#ifndef MATH_UTILS_H
#   define MATH_UTILS_H

#include <math.h>
#include "utils/macros.h"

#ifdef USE_NAMESPACE
namespace sploosh
{
#endif

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

inline int GCD(int a, int b)
{
    for(int c;b;c=a,a=b,b=c%b) ;
    return a;
}

inline int LCM(int a, int b) { return a/GCD(a,b)*b; }

template <typename T>
inline T M_CLAMP(T a, T minv, T maxv)
{
    return a <= minv ? minv : (a > maxv ? maxv : a);
}

template <typename T>
inline T M_GAUSSIAN(T v, T mu, T d) { return exp(-M_SQR((v - mu) / d)); }

#ifdef USE_NAMESPACE
}
#endif

#endif
