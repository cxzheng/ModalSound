#ifndef LINEARALGEBRA_EIG3_H
#   define LINEARALGEBRA_EIG3_H

/* Eigen-decomposition for symmetric 3x3 real matrices.
   Public domain, copied from the public domain Java library JAMA. */

/* 
 * Symmetric matrix A => eigenvectors in columns of V, corresponding
 * eigenvalues in d. 
 *
 * V[:][0] is the first eigenvector
 * V[:][1] is the second eigenvector
 * V[:][2] is the third eigenvector
 */
template<typename T>
void eigen_decomposition(const T A[3][3], T V[3][3], T d[3]);

template<>
void eigen_decomposition(const double A[3][3], double V[3][3], double d[3]);

#endif
