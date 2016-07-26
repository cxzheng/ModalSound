/*
 * =====================================================================================
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
 * 
 * -------------------------------------------------------------------------------------
 *
 *       Filename:  SymEigen.h
 *
 *    Description:  compute eigenvalues/eigenvectors of a symmetric matrix.
 *                  This class is a wrapper of MKL lapack routines
 *
 *        Version:  1.0
 *        Created:  07/22/12 22:48:54
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef SYM_EIGEN_INC
#   define SYM_EIGEN_INC

#include <vector>

template <typename T>
class SymEigen     
{
    public:
        SymEigen(int N);

        /* 
         * Solve the eigen problem 
         *
         * up: upper / lower part of the matrix is given
         * sm: compute the smallest n eigen-vectors. By default,
         *     it always computes the largest n eigen-vectors
         */
        int solve(T* A, int n, T* eval, T* evec, bool up = true, bool sm=false);

    private:
        int                 N_; /* size of the matrix */
        std::vector<T>      work_;
        std::vector<int>    iwork_;
        std::vector<int>    isuppz_;
        int                 lwork_;
        int                 liwork_;
};

template<>
SymEigen<double>::SymEigen(int N);

template<>
int SymEigen<double>::solve(double* A, int n, double* eval, double* evec, bool up, bool sm);

#endif
