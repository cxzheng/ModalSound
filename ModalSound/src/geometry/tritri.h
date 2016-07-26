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
 *       Filename:  tritri.h
 *
 *        Version:  1.0
 *        Created:  09/08/11 16:40:33
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef TRI_TRI_TEST_INC
#   define TRI_TRI_TEST_INC

//* function prototype *//
int tri_tri_overlap_test_3d(double* p);

int tri_tri_overlap_test_3d(const double p1[3], const double q1[3], 
                            const double r1[3], const double p2[3], 
                            const double q2[3], const double r2[3]);


int coplanar_tri_tri3d(const double p1[3], const double q1[3], 
                       const double r1[3], const double p2[3], 
                       const double q2[3], const double r2[3],
		       double N1[3], double N2[3]);


int tri_tri_overlap_test_2d(const double p1[2], const double q1[2], 
                            const double r1[2], const double p2[2],
                            const double q2[2], const double r2[2]);


/* coplanar returns whether the triangles are coplanar  
 * source and target are the endpoints of the segment of 
 * intersection if it exists) 
 */
int tri_tri_intersection_test_3d(const double p1[3], const double q1[3], const double r1[3], 
				 const double p2[3], const double q2[3], const double r2[3],
				 int * coplanar, 
				 double source[3],double target[3]);

#endif
