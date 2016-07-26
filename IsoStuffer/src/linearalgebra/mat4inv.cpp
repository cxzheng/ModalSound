/******************************************************************************
 *  File: mat4inv.cpp
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
#include "mat4inv.h"

void mat4d_invert(const double* mat, double* dst)
{
    double tmp[12];  /* temp array for pairs */
    double src[16];  /* array of transpose source matrix */
    double det;      /* determinant */
    
    /* transpose matrix */
    for (int i = 0; i < 4; i++) 
    {
        src[i] = mat[i*4];
        src[i + 4] = mat[i*4 + 1];
        src[i + 8] = mat[i*4 + 2];
        src[i + 12] = mat[i*4 + 3];
    }

    /* calculate pairs for first 8 elements (cofactors) */
    tmp[0] = src[10] * src[15];
    tmp[1] = src[11] * src[14];
    tmp[2] = src[9] * src[15];
    tmp[3] = src[11] * src[13];
    tmp[4] = src[9] * src[14];
    tmp[5] = src[10] * src[13];
    tmp[6] = src[8] * src[15];
    tmp[7] = src[11] * src[12];
    tmp[8] = src[8] * src[14];
    tmp[9] = src[10] * src[12];
    tmp[10] = src[8] * src[13];
    tmp[11] = src[9] * src[12];

    /* calculate first 8 elements (cofactors) */
    dst[0]  = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
    dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
    dst[1]  = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
    dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
    dst[2]  = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
    dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
    dst[3]  = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
    dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
    dst[4]  = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
    dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
    dst[5]  = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
    dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
    dst[6]  = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
    dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
    dst[7]  = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
    dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];

    /* calculate pairs for second 8 elements (cofactors) */
    tmp[0] = src[2]*src[7];
    tmp[1] = src[3]*src[6];
    tmp[2] = src[1]*src[7];
    tmp[3] = src[3]*src[5];
    tmp[4] = src[1]*src[6];
    tmp[5] = src[2]*src[5];
    tmp[6] = src[0]*src[7];
    tmp[7] = src[3]*src[4];
    tmp[8] = src[0]*src[6];
    tmp[9] = src[2]*src[4];
    tmp[10] = src[0]*src[5];
    tmp[11] = src[1]*src[4];

    /* calculate second 8 elements (cofactors) */
    dst[8]   = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
    dst[8]  -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
    dst[9]   = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
    dst[9]  -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
    dst[10]  = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
    dst[10] -= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
    dst[11]  = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
    dst[11] -= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
    dst[12]  = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
    dst[12] -= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
    dst[13]  = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
    dst[13] -= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
    dst[14]  = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
    dst[14] -= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
    dst[15]  = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
    dst[15] -= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];

    /* calculate determinant */
    det = src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];

    /* calculate matrix inverse */
    det = 1/det;
    for(int j = 0; j < 16; j++)
        dst[j] *= det;
}
