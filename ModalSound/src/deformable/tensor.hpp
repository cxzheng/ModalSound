#ifndef DEFORMABLE_TENSOR_INC
#   define DEFORMABLE_TENSOR_INC

#include "sc/Matrix3.hpp"

/*! 
 * Compute Right Cauchy-Green deformation tensor
 * C = F^T * F
 */
template <typename T>
inline void right_cauchy_green_deformation_tensor(
        const Matrix3<T>& F, Matrix3<T>& C)
{
    C = F.transpose() * F;
}

template <typename T>
inline void green_stain_tensor(const Matrix3<T>& C, Matrix3<T>& E)
{
    E = (C - Matrix3<T>::I) * 0.5;
}

#endif
