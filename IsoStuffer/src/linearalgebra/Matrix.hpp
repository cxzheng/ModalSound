/******************************************************************************
 *  File: Matrix.hpp
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
#ifndef LINEARALGEBRA_MATRIX_HPP
#   define LINEARALGEBRA_MATRIX_HPP

#include <assert.h>
#include <boost/multi_array.hpp>
#include <math.h>

template <typename T>
class Matrix : public boost::multi_array<T, 2>
{
    public:
        Matrix(int s, bool isSym = false):
                boost::multi_array<T, 2>(boost::extents[s][s], 
                                         boost::fortran_storage_order()), //: boost::c_storage_order()),
                m_isSymmetric(isSym)
        { }

        Matrix(int m, int n):
                boost::multi_array<T, 2>(boost::extents[m][n],
                                         boost::fortran_storage_order()), // : boost::c_storage_order()),
                m_isSymmetric(false)
        { }


        inline void zeros()
        {
            memset(this->data(), 0, sizeof(T)*this->num_elements());
        }

        inline void set(int m, int n, T v)
        {
            assert(m < this->shape()[0] && n < this->shape()[1]);
            (*this)[m][n] = v;

            if ( m_isSymmetric && m != n )
            {
                assert(m < this->shape()[1] && n < this->shape()[0]);
                (*this)[n][m] = v;
            }
        }

        bool check_symmetry() const
        {
            if ( this->shape()[0] != this->shape()[1] ) return false;

            for(int i = 0;i < this->shape()[0];++ i)
                for(int j = i+1;j < this->shape()[1];++ j)
                    if ( M_ABS((*this)[i][j] - (*this)[j][i]) > 1E-8 ) return false;
            return true;
        }

    private:
        bool                m_isSymmetric;
};

#endif
