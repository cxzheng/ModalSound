/******************************************************************************
 *  File: StdVectorOpt.hpp
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
#ifndef LINEARALGEBRA_STD_VECTOR_OPT_HPP
#   define LINEARALGEBRA_STD_VECTOR_OPT_HPP

#include <assert.h>
#include <vector>
#include <math.h>
#include "Vector3.hpp"

struct StdVectorOpt
{
    template <typename T>
    static inline void copy_from(std::vector<T>& dst, const std::vector<T>& src)
    {
        dst.resize(src.size());
        memcpy(&dst[0], &src[0], sizeof(T)*dst.size());
    }

    template <typename T>
    static void scale(std::vector<T>& vec, T s)
    {
        for(size_t i = 0;i < vec.size();++ i)
            vec[i] *= s;
    }

    template <typename T>
    static void scale(std::vector< Vector3<T> >& vec, T s)
    {
        for(size_t i = 0;i < vec.size();++ i)
            vec[i] *= s;
    }

    template <typename T>
    static void substract_from(std::vector<T>& out, const std::vector<T>& vec)
    {
        assert(out.size() == vec.size());
        for(size_t i = 0;i < out.size();++ i)
            out[i] -= vec[i];
    }

    template <typename T>
    static void add_to(std::vector<T>& out, const std::vector<T>& vec)
    {
        assert(out.size() == vec.size());
        for(size_t i = 0;i < out.size();++ i)
            out[i] += vec[i];
    }

    template <typename T>
    static void add(std::vector<T>& out, const std::vector<T>& a1, const std::vector<T>& a2)
    {
        assert(out.size() == a1.size() && out.size() == a2.size());
        for(size_t i = 0;i < out.size();++ i)
            out[i] = a1[i] + a2[i];
    }

    template <typename T>
    static T norm(const std::vector<T>& vec)
    {
        T ret = 0;
        for(size_t i = 0;i < vec.size();++ i)
            ret += vec[i]*vec[i];
        return sqrt(ret);
    }

    template <typename T>
    static T norm(const std::vector< Vector3<T> >& vec)
    {
        T ret = 0;
        for(size_t i = 0;i < vec.size();++ i)
            ret += vec[i].lengthSqr();
        return sqrt(ret);
    }

    template <typename T>
    static void axpy(std::vector< Vector3<T> >& out, T alpha, 
                     const std::vector< Vector3<T> >& vec)
    {
        assert(out.size() == vec.size());
        for(size_t i = 0;i < vec.size();++ i)
            out[i] += alpha * vec[i];
    }

};

#endif
