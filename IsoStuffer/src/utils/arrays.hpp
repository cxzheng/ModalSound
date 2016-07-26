/******************************************************************************
 *  File: arrays.hpp
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
 * arrays.hpp
 * author: Changxi Zheng (cxzheng@cs.cornell.edu)
 */
#ifndef UTILS_HPP
#   define UTILS_HPP

#include <boost/multi_array.hpp>
#include <string.h>
#include "math.hpp"

template<typename T, std::size_t D>
inline void zero_multi_array(boost::multi_array<T, D>& array)
{
    memset(array.data(), 0, sizeof(T)*array.num_elements());
}


template<typename T, std::size_t D>
struct MultiArrayOpt
{
    /*!
     * return the maximum norm of the multi-array
     */
    T sup_norm(const boost::multi_array<T, D>& array) const
    {
        const T* ptr = array.data();
        int      end = array.num_elements();
        T        ret = static_cast<T>(0);
        T        v;
    
        for(int i = 0;i < end;++ i,++ ptr)
            if ( (v = M_ABS(*ptr)) > ret ) ret = v;
        return ret;
    }

    /*!
     * TODO: Use blas to accelerate it
     *       it can also be parallized
     */
    T dot(const boost::multi_array<T, D>& a,
            const boost::multi_array<T, D>& b) const
    {
        T ret = 0;
        const T* ptrA = a.data();
        const T* ptrB = b.data();
    
        for(size_t i = 0;i < a.num_elements();++ i)
            ret += ptrA[i]*ptrB[i];
        return ret;
    }
    
    /*!
     * <aout> = <aout> + scale*<ain>
     */
    void scaleAdd(boost::multi_array<T, D>& aout,
            const boost::multi_array<T, D>& ain, const T scale) const
    {
        const T* ptrIn = ain.data();
        T* ptrOut = aout.data();
        for(size_t i = 0;i < ain.num_elements();++ i)
            ptrOut[i] += scale*ptrIn[i];
    }

    /*!
     * <aout> = scale*<aout> + <ain>
     */
    void selfScaleAdd(boost::multi_array<T, D>& aout,
            const T scale, const boost::multi_array<T, D>& ain) const
    {
        const T* ptrIn = ain.data();
        T* ptrOut = aout.data();
        for(size_t i = 0;i < ain.num_elements();++ i)
            ptrOut[i] = ptrOut[i]*scale + ptrIn[i];
    }
};

#ifdef USE_OPENMP
/*!
 * use OpenMP for parallization 
 */
template<typename T, std::size_t D>
struct ParallelMultiArrayOpt
{
    /*!
     * return the maximum norm of the multi-array
     */
    T sup_norm(const boost::multi_array<T, D>& array) const
    {
        const T* ptr = array.data();
        int      end = array.num_elements();
        T        ret = static_cast<T>(0);
        T        v;
    
        for(int i = 0;i < end;++ i,++ ptr)
            if ( (v = M_ABS(*ptr)) > ret ) ret = v;
        return ret;
    }

    T dot(const boost::multi_array<T, D>& a,
            const boost::multi_array<T, D>& b) const
    {
        T ret = 0;
        const T* ptrA = a.data();
        const T* ptrB = b.data();
        int end = a.num_elements();
    
        #pragma omp parallel for reduction(+:ret) schedule(dynamic, 8000) \
                default(none) shared(end, ptrA, ptrB)
        for(int i = 0;i < end;++ i)
            ret += ptrA[i]*ptrB[i];
        return ret;
    }
    
    /*!
     * <aout> = <aout> + scale*<ain>
     */
    void scaleAdd(boost::multi_array<T, D>& aout,
            const boost::multi_array<T, D>& ain, const T scale) const
    {
        const T* ptrIn = ain.data();
        T* ptrOut = aout.data();
        #pragma omp parallel for default(none) schedule(dynamic, 8000) \
                shared(ain, ptrOut, ptrIn)
        for(int i = 0;i < (int)ain.num_elements();++ i)
            ptrOut[i] += scale*ptrIn[i];
    }

    /*!
     * <aout> = scale*<aout> + <ain>
     */
    void selfScaleAdd(boost::multi_array<T, D>& aout,
            const T scale, const boost::multi_array<T, D>& ain) const
    {
        const T* ptrIn = ain.data();
        T* ptrOut = aout.data();
        #pragma omp parallel for default(none) schedule(dynamic, 8000) \
                shared(ain, ptrOut, ptrIn)
        for(int i = 0;i < (int)ain.num_elements();++ i)
            ptrOut[i] = ptrOut[i]*scale + ptrIn[i];
    }
};
#endif

#endif
