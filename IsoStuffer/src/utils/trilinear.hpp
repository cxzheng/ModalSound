/******************************************************************************
 *  File: trilinear.hpp
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
#ifndef UTILS_TRILINEAR_HPP
#   define UTILS_TRILINEAR_HPP

/*
 * Trilinear interpolation:
 * val[2][2][2] stores the value at eight adjacent (voxel) points.
 * val[Z][Y][X] 
 */
template <typename T, typename V>
inline T trilinear_interpolate(const T val[][2][2], const V* alpha)
{
    return ((val[0][0][0]*(1. - alpha[0]) + val[0][0][1]*alpha[0]) *
            (1. - alpha[1]) +
            (val[0][1][0]*(1. - alpha[0]) + val[0][1][1]*alpha[0]) *
            alpha[1]) * (1. - alpha[2]) + 
           ((val[1][0][0]*(1. - alpha[0]) + val[1][0][1]*alpha[0]) *
            (1. - alpha[1]) + 
            (val[1][1][0]*(1. - alpha[0]) + val[1][1][1]*alpha[0]) *
            alpha[1]) * alpha[2];
}

#endif
