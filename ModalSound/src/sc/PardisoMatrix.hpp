/******************************************************************************
 *  File: PardisoMatrix.hpp
 *  A sparse matrix used for Intel MKL PARDISO and DSS solver
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
#ifndef SC_PARDISO_MATRIX_HPP
#   define SC_PARDISO_MATRIX_HPP

#include <string.h>
#include <map>
#include <assert.h>
#include "DiagonalMatrix.hpp"
#include "Matrix.hpp"
#include "Vector3.hpp"

/*!
 * Sparse Matrix
 * The data structure of this sparse matrix is particularly optimized
 * for intel MKL PARDISO or DSS solver call. 
 *
 * NOTE: it assumes the diagonal elements is always non-zero
 */

struct PardisoMatrixIO;

template <typename T>
class PardisoMatrix
{
    friend struct PardisoMatrixIO;

    public:
        PardisoMatrix(int m, int n):
                numCols_(n), numRows_(m), isSymm_(false),
                isPosDef_(false)
        {  rows_.resize(m); }

        PardisoMatrix(int s, bool isSym = false):
                numCols_(s), numRows_(s), isSymm_(isSym),
                isPosDef_(false)
        {  rows_.resize(s); }

        PardisoMatrix(int s, bool isSym, bool isPosDef):
                numCols_(s), numRows_(s), isSymm_(isSym),
                isPosDef_(isPosDef)
        {  rows_.resize(s); }

        PardisoMatrix():numCols_(0), numRows_(0), isSymm_(false),
                isPosDef_(false)
        {}

        void clear()
        {
            for(size_t i = 0;i < rows_.size();++ i) rows_[i].clear();
        }

        void resize(int s, bool isSym = false)
        {
            numCols_ = numRows_ = s;
            isSymm_ = isSym;
            rows_.resize(s);
            clear();
        }

        void resize(int s, bool isSym, bool isPosDef)
        {
            numCols_ = numRows_ = s;
            isSymm_ = isSym;
            isPosDef_ = isPosDef;
            rows_.resize(s);
            clear();
        }

        void add(const PardisoMatrix<T>& mat)
        {
            if ( data_.size() != mat.data_.size() ||
                 rowIdx_.size() != mat.rowIdx_.size() )
            {
                fprintf(stderr, "ERROR: Two Pardiso Matrices for scaleAdd call have different sparsity pattern\n");
                exit(1);
            }

            for(size_t i = 0;i < data_.size();++ i)
                data_[i] += mat.data_[i];
        }

        /*!
         * Axpy operator: compute this = this + beta * mat
         * It assumes the this matrix and the mat have the same sparsity pattern.
         */
        void axpy(T beta, const PardisoMatrix<T>& mat)
        {
            if ( data_.size() != mat.data_.size() ||
                 rowIdx_.size() != mat.rowIdx_.size() )
            {
                fprintf(stderr, "ERROR: Two Pardiso Matrices for scaleAdd call have different sparsity pattern\n");
                exit(1);
            }

            for(size_t i = 0;i < data_.size();++ i)
                data_[i] += beta * mat.data_[i];
        }

        void axpy(T beta, const DiagonalMatrix<T>& mat)
        {
            if ( mat.size() != numCols_ || numCols_ != numRows_ )
            {
                fprintf(stderr, "ERROR: the size of diagonal matrix is inconsistent with this PardisoMatrix\n");
                exit(1);
            }

#ifdef USE_OPENMP
#warning ------------ OPENMP is enabled -----------
            #pragma omp parallel for default(none) schedule(dynamic, 8000) shared(beta, mat)
#endif
            for(int i = 0;i < numCols_;++ i)
                *(rows_[i][i]) += beta * mat[i];
        }

        void axpy(T beta, const DiagonalMatrix< Vector3<T> >& mat)
        {
            if ( mat.size()*3 != numCols_ || numCols_ != numRows_ )
            {
                fprintf(stderr, "ERROR: the size of diagonal matrix is inconsistent with this PardisoMatrix\n");
                exit(1);
            }

            const T* ptr = (const T*)mat.data();
#ifdef USE_OPENMP
            #pragma omp parallel for default(none) schedule(dynamic, 8000) shared(beta, ptr)
#endif
            for(int i = 0;i < numCols_;++ i)
                *(rows_[i][i]) += beta * ptr[i];
        }

        /*! 
         * indicate the mat[nr][nc] is nonzero
         * the benifit is, for example,
         * as long as we don't remesh, and the mesh isn't split during
         * the simulation, the sparity pattern of stiffness matrix 
         * doesn't change.
         */
        void set_nonzero(int nr, int nc);
        void generate_pattern();        // $$TESTED

        void add(int nr, int nc, T v)
        {
            if ( isSymm_ && nr > nc )
            {
                /*
                 * Here we assume mat[nr][nc] += v only applies to the
                 * upper trangle part of the matrix when this is a 
                 * symmetric matrix
                 */
                return;
            }
            assert(rows_[nr].count(nc));
            *(rows_[nr][nc]) += v;
        }

        /*! make all the elements be zero */
        void zeros()
        {
            memset(&data_[0], 0, sizeof(T)*data_.size());
        }

        void multiply(const std::vector< Vector3<T> >& in, 
                      std::vector< Vector3<T> >& out) const;

        void multiply(const T* in, T* out) const;

        bool check_symmetry();

        void dump() const;
        void to_matrix(Matrix<T>& mat) const;

        /* ================ Retrival Methods =============== */
        int num_rows() const
        { return numRows_; }
        int num_cols() const
        { return numCols_; }
        int num_nonzeros() const
        { return data_.size(); }

        const T* data() const
        { return &data_[0]; }

        const int* col_indices() const
        { return &cols_[0]; }

        const int* row_indices() const
        { return &rowIdx_[0]; }

        T* data() 
        { return &data_[0]; }

        int* col_indices() 
        { return &cols_[0]; }

        int* row_indices()
        { return &rowIdx_[0]; }

        bool is_symmetric() const 
        { return isSymm_; }

        bool is_positive_definite() const
        { return isPosDef_; }

        bool& is_positive_definite() 
        { return isPosDef_; }

        bool operator == (const PardisoMatrix<T>& m2) const;

    private:
        typedef std::map<int, T*> RowMap;

        int                 numCols_;
        int                 numRows_;
        int                 numNonZero_;
        std::vector<RowMap> rows_;
        /*! indicate if it is symmetric matrix */
        bool                isSymm_;
        bool                isPosDef_;

        /* the data structure for Pardiso storage format */
        std::vector<T>      data_;
        std::vector<int>    cols_;
        std::vector<int>    rowIdx_;
};

///////////////////////////////////////////////////////////////////////////////
template <typename T>
bool PardisoMatrix<T>::check_symmetry()
{
    if ( numCols_ != numRows_ ) return false;
    if ( isSymm_ ) return true;

    for(int i = 0;i < numRows_;++ i)
    {
        typename RowMap::iterator end = rows_[i].end();
        for(typename RowMap::iterator it = rows_[i].begin();
                it != end;++ it)
        {
            if ( it->first < i ) continue;
            if ( !rows_[it->first].count(i) || 
                 M_ABS(*rows_[it->first][i] - *it->second) > 
                 1E-10*M_MAX(M_ABS(*rows_[it->first][i]), M_ABS(*it->second)) )
                return false;
        }
    }
    return true;
}

/*
 * If the matrix is symmetric, it stores only the upper part of it
 */
template <typename T>
void PardisoMatrix<T>::set_nonzero(int nr, int nc)
{
    if ( isSymm_ && nr > nc ) 
        rows_[nc][nr] = NULL;
    else
        rows_[nr][nc] = NULL;
}

/*
 * NOTE that in order to ease to work with Fortran routines, the indices
 *      in both cols_ and rowIdx_ are all 1-based.
 */
template <typename T>
void PardisoMatrix<T>::generate_pattern()
{
    data_.clear();
    cols_.clear();
    rowIdx_.clear();

    // # of nonzero elements
    numNonZero_ = 0;
    for(int i = 0;i < numRows_;++ i)
        numNonZero_ += rows_[i].size();
    data_.resize(numNonZero_);

    for(int i = 0, ptr = 0;i < numRows_;++ i)
    {
        if ( rows_.empty() )
        {
            fprintf(stderr, "Warning: With zeros in the %d-th row\n, cannot generate sparsity pattern for PARDISO call", i);
            return;
        }

        typename RowMap::iterator it  = rows_[i].begin();
        typename RowMap::iterator end = rows_[i].end();
        rowIdx_.push_back(ptr + 1);
        for(;it != end;++ it)
        {
            cols_.push_back(it->first + 1);    // 1-based column index
            it->second = &(data_[ptr]);
            ++ ptr;
        }
    }

    rowIdx_.push_back(data_.size() + 1);      // 1-based row index
}

template <typename T>
void PardisoMatrix<T>::multiply(const std::vector< Vector3<T> >& in, 
                                std::vector< Vector3<T> >& out) const
{
    assert(in.size()*3 == numCols_);
    if ( out.size() != numRows_ ) out.resize(numRows_);

    const T* inPtr = (const T*)&in[0];
    T* outPtr = (T*)out[0];
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic) shared(inPtr, outPtr)
#endif
    for(int i = 0;i < numRows_;++ i)
    {
        outPtr[i] = 0;
        typename RowMap::const_iterator end = rows_[i].end();
        for(typename RowMap::const_iterator it = rows_[i].begin();
                it != end;++ it)
            outPtr[i] += inPtr[it->first]*(*(it->second));
    }
}

template <typename T>
void PardisoMatrix<T>::multiply(const T* in, T* out) const
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic) shared(inPtr, outPtr)
#endif
    for(int i = 0;i < numRows_;++ i)
    {
        out[i] = 0;
        typename RowMap::const_iterator end = rows_[i].end();
        for(typename RowMap::const_iterator it = rows_[i].begin();
                it != end;++ it)
            out[i] += in[it->first]*(*(it->second));
    }
}

template <typename T>
void PardisoMatrix<T>::to_matrix(Matrix<T>& mat) const
{
    mat.zeros();

    int rowId = 0;
    for(size_t i = 0;i < data_.size();++ i)
    {
        int colId = cols_[i] - 1;
        if ( rowIdx_[rowId+1] == i+1 ) ++ rowId;
        mat.set(colId, rowId, data_[i]);
    }
}

template <typename T>
void PardisoMatrix<T>::dump() const
{
    using namespace std;

    cout << "data: { ";
    for(int i = 0;i < data_.size();++ i)
        cout << data_[i] << ' ';
    cout << '}' << endl;
    cout << "col: { ";
    for(int i = 0;i < cols_.size();++ i)
        cout << cols_[i] << ' ';
    cout << '}' << endl;
    cout << "row: { ";
    for(int i = 0;i < rowIdx_.size();++ i)
        cout << rowIdx_[i] << ' ';
    cout << '}' << endl;
}

template <>
void PardisoMatrix<double>::add(const PardisoMatrix<double>& mat);

template<>
void PardisoMatrix<double>::axpy(double beta, const PardisoMatrix<double>& mat);

template<>
bool PardisoMatrix<double>::operator == (const PardisoMatrix<double>& m2) const;

template<>
bool PardisoMatrix<float>::operator == (const PardisoMatrix<float>& m2) const;

#endif
