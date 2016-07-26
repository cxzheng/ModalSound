#ifndef HELMHOLTZ_BASIS_HPP
#   define HELMHOLTZ_BASIS_HPP

#include <complex>
#include <gsl/gsl_sf.h>

#include "geometry/Point3.hpp"
#include "utils/math.hpp"

namespace HelmholtzBasis
{

template <typename T>
inline std::complex<T> spherical_hankel_1st(int n, T kr)
{
    return std::complex<T>(gsl_sf_bessel_jl(n, kr), gsl_sf_bessel_yl(n, kr));
}

template <typename T>
inline std::complex<T> spherical_harmonics(int m, int n, T theta, T phi)
{
    return m & 1 ? -gsl_sf_legendre_sphPlm(n, M_ABS(m), cos(theta)) * std::complex<T>(cos(m*phi), sin(m*phi)) :
                    gsl_sf_legendre_sphPlm(n, M_ABS(m), cos(theta)) * std::complex<T>(cos(m*phi), sin(m*phi));
}

template <typename T>
static inline T Amn(int m, int n)
{
    return n < M_ABS(m) ? 0 : (T)(sqrt(double((m+n+1)*(n-m+1)) / double((2*n+1)*(2*n+3))));
}

template <typename T>
static inline T Bmn(int m, int n)
{
    return n < M_ABS(m) ? 0 :
        (m >= 0 ? (T)(sqrt(double((n-m-1)*(n-m)) / double((2*n-1)*(2*n+1)))) : 
                 -(T)(sqrt(double((n-m-1)*(n-m)) / double((2*n-1)*(2*n+1)))));
}

template <typename T>
static std::complex<T> regular_basis(int m, int n, T K, const Point3<T>& sc)
{
    if ( n < M_ABS(m) ) return std::complex<T>(0, 0);
    return gsl_sf_bessel_jl(n, K*sc.r) * spherical_harmonics(m, n, sc.theta, sc.phi);
}

template <typename T>
static std::complex<T> regular_basis_dir_deriv(int m, int n, T K, 
        const Point3<T>& sc, const Vector3<T>& dir)
{
    using namespace std;
    const complex<T> A = Bmn<T>(-m-1, n+1) * regular_basis(m+1, n+1, K, sc) - 
                         Bmn<T>(m,    n)   * regular_basis(m+1, n-1, K, sc);
    const complex<T> B = Bmn<T>(m-1, n+1)  * regular_basis(m-1, n+1, K, sc) -
                         Bmn<T>(-m,  n)    * regular_basis(m-1, n-1, K, sc);
    complex<T> rx = (0.5*A + 0.5*B);
    complex<T> ry = complex<T>((A.imag() - B.imag()) * 0.5,
                               (B.real() - A.real()) * 0.5);
    complex<T> rz = Amn<T>(m, n-1) * regular_basis(m, n-1, K, sc) -
                    Amn<T>(m, n)   * regular_basis(m, n+1, K, sc);
    return (rx*dir.x + ry*dir.y + rz*dir.z)*K;
}

template <typename T>
static std::complex<T> singular_basis(int m, int n, T K, const Point3<T>& sc)
{
    if ( n < M_ABS(m) ) return std::complex<T>(0, 0);
    return spherical_hankel_1st(n, K*sc.r) * spherical_harmonics(m, n, sc.theta, sc.phi);
}

//// ------ to compute the derivative of basis function ------
template <typename T>
static T basis_deri_A(int m, int n)
{
    if ( n < M_ABS(m) ) return (T)0;
    int k = n+1;
    T numer = (T)((k-m)*(k+m));
    k = 2*n + 1;
    T denom = (T)(k*(k+2));     // (2n+1)(2n+3)
    return sqrt(numer / denom);
}

template <typename T>
static T basis_deri_B(int m, int n)
{
    if ( n < M_ABS(m) || m == n ) return (T)0;
    int k = n-m;
    T numer = (T)((k-1)*k);     // (n-m-1)(n-m)
    k = 2*n - 1;
    T denom = (T)(k*(k+2));     // (2n-1)(2n+1)
    T ret = sqrt(numer / denom);
    return m >= 0 ? ret : -ret;
}

}

#endif
