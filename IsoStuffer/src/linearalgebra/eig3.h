/******************************************************************************
 *  File: eig3.h
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
void eigen_decomposition(const double A[3][3], double V[3][3], double d[3]);

#endif
