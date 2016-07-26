#include "PardisoMatrix.hpp"
#include <stdio.h>

#include <map>
#ifdef USE_MKL
#   include <mkl.h>
#else
#   if defined(__APPLE__) || defined(MACOSX)
#       include <vecLib/cblas.h>
#   else
#       include <cblas.h>
#   endif
#endif

#include "utils/math.hpp"

template <>
void PardisoMatrix<double>::add(const PardisoMatrix<double>& mat)
{
    if ( data_.size() != mat.data_.size() ||
         rowIdx_.size() != mat.rowIdx_.size() )
    {
        fprintf(stderr, "ERROR: Two Pardiso Matrices for scaleAdd call have different sparsity pattern\n");
        exit(1);
    }

    cblas_daxpy(data_.size(), 1, &mat.data_[0], 1, &data_[0], 1);
}

template<>
void PardisoMatrix<double>::axpy(double beta, const PardisoMatrix<double>& mat)
{
    if ( data_.size() != mat.data_.size() ||
         rowIdx_.size() != mat.rowIdx_.size() )
    {
        fprintf(stderr, "ERROR: Two Pardiso Matrices for scaleAdd call have different sparsity pattern\n");
        exit(1);
    }

    cblas_daxpy(data_.size(), beta, &mat.data_[0], 1, &data_[0], 1);
}

template<>
bool PardisoMatrix<double>::operator == (const PardisoMatrix<double>& m2) const
{
    if ( cols_.size() != m2.cols_.size() ||
         rowIdx_.size() != m2.rowIdx_.size() ||
         data_.size() != m2.data_.size() )
        return false;

    for(int i = (int)cols_.size()-1;i >= 0;-- i)
        if ( cols_[i] != m2.cols_[i] ) return false;
    for(int i = (int)rowIdx_.size()-1;i >= 0;-- i)
        if ( rowIdx_[i] != m2.rowIdx_[i] ) return false;
    for(int i = (int)data_.size()-1;i >= 0;-- i)
        if ( M_ABS(data_[i]-m2.data_[i]) > 1E-12 ) return false;
    return true;
}

template<>
bool PardisoMatrix<float>::operator == (const PardisoMatrix<float>& m2) const
{
    if ( cols_.size() != m2.cols_.size() ||
         rowIdx_.size() != m2.rowIdx_.size() ||
         data_.size() != m2.data_.size() )
        return false;

    for(int i = (int)cols_.size()-1;i >= 0;-- i)
        if ( cols_[i] != m2.cols_[i] ) return false;
    for(int i = (int)rowIdx_.size()-1;i >= 0;-- i)
        if ( rowIdx_[i] != m2.rowIdx_[i] ) return false;
    for(int i = (int)data_.size()-1;i >= 0;-- i)
        if ( M_ABS(data_[i]-m2.data_[i]) > (float)1E-8 ) return false;
    return true;
}

