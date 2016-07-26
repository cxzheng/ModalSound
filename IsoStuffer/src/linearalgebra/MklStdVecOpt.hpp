/******************************************************************************
 *  File: MklStdVecOpt.hpp
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
#ifndef LINEARALGEBRA_MKL_STD_VEC_OPT_HPP
#   define LINEARALGEBRA_MKL_STD_VEC_OPT_HPP

#include <mkl.h>
#include <assert.h>
#include <vector>
#include <math.h>
#include "Vector3.hpp"

struct MklArrayOpt
{
    static inline void add_to(double* out, const double* vec, int len)
    {
        cblas_daxpy(len, 1, vec, 1, out, 1);
    }

    static inline void copy_from(double* dst, const double* src, int len)
    {
        cblas_dcopy(len, src, 1, dst, 1);
    }

    static inline void scale(double* vec, double s, int len)
    {
        cblas_dscal(len, s, vec, 1);
    }

    static inline void substract_from(double* out, const double* vec, int len)
    {
        cblas_daxpy(len, -1, vec, 1, out, 1);
    }

    /*
     * out = out + alpha*vec
     */
    static inline void axpy(double* out, double alpha, const double* vec, int len)
    {
        cblas_daxpy(len, alpha, vec, 1, out, 1);
    }

    static inline double norm(const double* vec, int len)
    {
        return cblas_dnrm2(len, vec, 1);
    }

};

struct MklStdVecOpt
{
    template <typename T>
    static inline void copy_from(std::vector<T>& dst, const std::vector<T>& src)
    {
        dst.resize(src.size());
        memcpy(&dst[0], &src[0], sizeof(T)*dst.size());
    }
    
    static inline void copy_from(std::vector<Vector3d>& dst, 
                                 const std::vector<Point3d>& src)
    {
        dst.resize(src.size());
        cblas_dcopy((int)src.size()*3, (const double*)&src[0], 1, (double*)&dst[0], 1);
    }

    static inline void scale(std::vector<double>& vec, double s)
    {
        cblas_dscal(vec.size(), s, &vec[0], 1);
    }

    static inline void scale(std::vector< Vector3<double> >& vec, double s)
    {
        cblas_dscal(vec.size()*3, s, (double *)&vec[0], 1);
    }

    static inline void substract_from(std::vector< Point3<double> >& out, 
                                      const std::vector< Vector3<double> >& vec)
    {
        assert(out.size() == vec.size());
        cblas_daxpy(out.size()*3, -1, (double *)&vec[0], 1, (double *)&out[0], 1);
    }

    static inline void substract_from(std::vector< Vector3<double> >& out, 
                                      const std::vector< Vector3<double> >& vec)
    {
        assert(out.size() == vec.size());
        cblas_daxpy(out.size()*3, -1, (double *)&vec[0], 1, (double *)&out[0], 1);
    }

    static inline void add_to(std::vector< Vector3<double> >& out, 
                              const std::vector< Vector3<double> >& vec)
    {
        assert(out.size() == vec.size());
        cblas_daxpy(out.size()*3, 1, (double *)&vec[0], 1, (double *)&out[0], 1);
    }

    static inline double norm(const std::vector< Vector3<double> >& vec)
    {
        return cblas_dnrm2(vec.size()*3, (const double *)&vec[0], 1);
    }

    static inline void axpy(std::vector< Vector3<double> >& out, double alpha, 
                     const std::vector< Vector3<double> >& vec)
    {
        assert(out.size() == vec.size());
        cblas_daxpy(out.size()*3, alpha, (const double *)&vec[0], 1, (double *)&out[0], 1);
    }
    
};

#endif
