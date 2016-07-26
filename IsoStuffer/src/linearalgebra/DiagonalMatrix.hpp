/******************************************************************************
 *  File: DiagonalMatrix.hpp
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
#ifndef LINEARALGEBRA_DIAGONAL_MATRIX_HPP
#   define LINEARALGEBRA_DIAGONAL_MATRIX_HPP

#include <vector>

/*!
 * Diagonal matrix. i.e. The mass matrix in unreduced model
 */
template <typename T>
class DiagonalMatrix
{
    public:
        int size() const { return m_data.size(); }

        DiagonalMatrix<T>& operator = (const std::vector<T>& vec)
        {
            m_data = vec;
            return *this;
        }

        void resize(int n) 
        {
            m_data.resize(n); 
        }

        T operator [] (int n) const
        {
            assert(n >= 0 && n < (int)m_data.size());
            return m_data[n];
        }

        T& operator [] (int n) 
        {
            assert(n >= 0 && n < (int)m_data.size());
            return m_data[n];
        }

        void multiply(const std::vector<T>& in, std::vector<T>& out) const
        {
            assert(in.size() == m_data.size());
            if ( out.size() != in.size() ) out.resize(in.size());

#ifdef USE_OPENMP
            #pragma omp parallel for default(none) schedule(dynamic, 5000) shared(out, in)
#endif
            for(int i = 0;i < (int)m_data.size();++ i)
                out[i] = in[i] * m_data[i];
        }

        void multiply(const T* in, T* out) const
        {
#ifdef USE_OPENMP
            #pragma omp parallel for default(none) schedule(dynamic, 5000) shared(out, in)
#endif
            for(int i = 0;i < (int)m_data.size();++ i)
                out[i] = in[i] * m_data[i];
        }

        const T* data() const
        { return &m_data[0]; }

    private:
        std::vector<T>  m_data;
};

#endif
