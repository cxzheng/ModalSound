#ifndef SPHERICAL_HANKEL_HPP
#   define SPHERICAL_HANKEL_HPP

#include <complex>
#include "geometry/Point3.hpp"
#include "Vector3.hpp"
#include "generic/precision_type.hpp"

/*!
 * out.x ==> r
 * out.y ==> theta
 * out.z ==> phi
 */
template<typename T>                /* $TESTED$ */
inline void cartesian_to_spherical(const Point3<T>& from,
        const Point3<T>& to, Point3<T>& out)
{
    Vector3<T> v = to - from;
    out.r = v.length();

    if ( out.r < PrecisionType<T>::EPS )
    {
        out.zero();
        fprintf(stderr, "ERROR: <cartesian_to_spherical> r is too small! %lf\n", out.r);
        return;
    }

    out.theta = acos(v.z / out.r); 
    out.phi   = atan2(v.y, v.x);
}

/*
 * \return v.x^2 + v.y^2, this value is useful when computing the 
 *         gradient/directional derivative.
 */
template<typename T>                /* $TESTED$ */
inline T cartesian_to_spherical(const Vector3<T>& v, Point3<T>& out)
{
    out.r = v.length();

    if ( out.r < PrecisionType<T>::EPS )
    {
        out.zero();
        fprintf(stderr, "ERROR: <cartesian_to_spherical> r is too small! %lf\n", out.r);
        return 0;
    }

    T ret;
    if ( (ret = v.x*v.x + v.y*v.y) < PrecisionType<T>::EPS )
        fprintf(stderr, "ERROR: <cartesian_to_spherical> v.x^2+v.y^2 is too small! %lf\n", ret);

    out.theta = acos(v.z / out.r); 
    out.phi   = atan2(v.y, v.x);

    return ret;
}

#endif
