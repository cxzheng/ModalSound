#ifndef SC_MATRIX_HPP
#   define SC_MATRIX_HPP

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
                isSymmetric_(isSym)
        { }

        Matrix(int m, int n):
                boost::multi_array<T, 2>(boost::extents[m][n],
                                         boost::fortran_storage_order()), // : boost::c_storage_order()),
                isSymmetric_(false)
        { }


        inline void zeros()
        {
            memset(this->data(), 0, sizeof(T)*this->num_elements());
        }

        inline void set(int m, int n, T v)
        {
            assert(m < this->shape()[0] && n < this->shape()[1]);
            (*this)[m][n] = v;

            if ( isSymmetric_ && m != n )
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
        bool                isSymmetric_;
};

#endif
