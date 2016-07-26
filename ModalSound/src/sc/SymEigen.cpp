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
 *       Filename:  SymEigen.cpp
 *
 *        Version:  1.0
 *        Created:  07/22/12 23:40:54
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#include "SymEigen.hpp"
#include <mkl_lapack.h>
#include "utils/macros.h"
#include "utils/term_msg.h"

template<>
SymEigen<double>::SymEigen(int N):N_(N), lwork_(-1), liwork_(-1)
{
    isuppz_.resize(2*N);
    work_.resize(16);
    iwork_.resize(16);

    char jobz  = 'V';
    char range = 'I';
    char uplo  = 'U';
    int  il    = 1;
    int  iu    = N;
    double abstol = 0.;
    int  m, info;

    dsyevr(&jobz, &range, &uplo, &N, NULL, &N, NULL, NULL, &il, &iu,
           &abstol, &m, NULL, NULL, &N, &(isuppz_[0]), &(work_[0]),
           &lwork_, &(iwork_[0]), &liwork_, &info);
    if ( info ) 
    {
        PRINT_ERROR("Error: dsyevr has failed. info=%d\n", info);
        SHOULD_NEVER_HAPPEN(1);
    }

    liwork_ = iwork_[0];
    lwork_  = static_cast<int>(work_[0]);
    iwork_.resize(liwork_);
    work_.resize(lwork_);
    PRINT_MSG("Allocated %d int and %d double for workspace\n", liwork_, lwork_);
}

template<>
int SymEigen<double>::solve(double* A, int n, double* eval, double* evec, bool up, bool sm)
{
    char jobz = 'V';
    char range = 'I';
    char uplo = up ? 'U' : 'L';
    double abstol = 0;
    int il, iu, m, info;
    if ( unlikely(sm) )
    {
        il = 1; iu = n;
    }
    else
    {
        il = N_-n+1; iu = N_;
    }

    dsyevr(&jobz, &range, &uplo, &N_, A, &N_, NULL, NULL, &il, &iu,
           &abstol, &m, eval, evec, &N_, &(isuppz_[0]), &(work_[0]),
           &lwork_, &(iwork_[0]), &liwork_, &info);
    if ( info )
    {
        PRINT_ERROR("Error: dsyevr has failed. info=%d\n", info);
        return info;
    }
    return 0;
}

