/******************************************************************************
 *  File: DirectSparseSolver.hpp
 *  A linear solver using Intel MKL DSS routines
 *  Copyright (c) 2007 by Changxi Zheng
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
#ifndef SC_DIRECT_SPARSE_SOLVER_INC
#   define SC_DIRECT_SPARSE_SOLVER_INC

#include <mkl_dss.h>
#include "PardisoMatrix.hpp"

/*!
 * This is basically a wrap of MKL's Direct Sparse Solver
 *
 * It is used for the case where A is the same, and there are multiple b for
 * the linear system Ax = b. Therefore we can precompute the factorization of
 * A.
 */
class DirectSparseSolver
{
    public:
        DirectSparseSolver():opt_(MKL_DSS_DEFAULTS), allocated_(false)
        { }

        ~DirectSparseSolver()
        {
            if ( allocated_ ) dss_delete(handle_, opt_);
        }

        /*!
         * This method first free any previous allocated memory
         * - call dss_create to create a new handle for solver
         * - call dss_define_structure to tell the matrix structure
         *   NOTE: only MKL_DSS_SYMMETRIC / MKL_DSS_NON_SYMMETRIC would be passed
         *         to define the matrix structure
         *
         * - call dss_reorder to compute permutation matrix
         * - call dss_factor_real to factorize the given matrix
         *   NOTE: only MKL_DSS_POSITIVE_DEFINITE / MKL_DDS_INDEFINITE would be
         *         passed for factoring
         */
        void load_matrix(const PardisoMatrix<double>* A);

        /* 
         * solve the system with pre-loaded matrix A
         * it solves Ax = b, and the results is in x
         */
        void solve(const double* b, double *x);

        int& options()
        {  return opt_; }
        int options() const
        {  return opt_; }

    private:
        _MKL_DSS_HANDLE_t   handle_;
        int                 opt_;
        bool                allocated_;
};
#endif
