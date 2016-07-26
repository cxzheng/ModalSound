/******************************************************************************
 *  File: CubicMeshGen.h
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
#ifndef CUBIC_MESH_GEN_HPP
#   define CUBIC_MESH_GEN_HPP

#include "constants.h"
#include "utils/macros.h"
#include "geometry/FixVtxTetMesh.hpp"

struct CubicMeshGen
{
    static REAL dx, dy, dz;
    static int NX, NY, NZ;
    static REAL SIZE;

    static int load_mesh(const char*, FixVtxTetMesh<REAL>& mesh)
    {
        mesh.clear();
        for(int i = 0;i <= NY;++ i)
        for(int j = 0;j <= NZ;++ j)
        for(int k = 0;k <= NX;++ k)
            mesh.add_vertex(Point3<REAL>(k*SIZE + dx, i*SIZE + dy, j*SIZE + dz));

        mesh.set_fixed_vtx((NX+1)*(NZ+1));
        for(int i = 0;i < NY;++ i)
        for(int j = 0;j < NZ;++ j)
        for(int k = 0;k < NX;++ k)
        {
            if ( (i+j+k) % 2 )
            {
                mesh.add_tet(get_id(k+1,i+1,j+1), get_id(k+1,i+1,j), 
                             get_id(k+1,i,  j+1), get_id(k  ,i+1,j+1));
                mesh.add_tet(get_id(k  ,i+1,j  ), get_id(k+1,i+1,j),
                             get_id(k  ,i+1,j+1), get_id(k  ,i  ,j));
                mesh.add_tet(get_id(k+1,i  ,j  ), get_id(k+1,i+1,j),
                             get_id(k  ,i  ,j  ), get_id(k+1,i  ,j+1));
                mesh.add_tet(get_id(k  ,i  ,j+1), get_id(k  ,i+1,j+1),
                             get_id(k+1,i  ,j+1), get_id(k  ,i  ,j));
                mesh.add_tet(get_id(k  ,i+1,j+1), get_id(k+1,i+1,j),
                             get_id(k+1,i  ,j+1), get_id(k  ,i  ,j));
            }
            else
            {
                mesh.add_tet(get_id(k+1,i+1,j  ), get_id(k+1,i  ,j),
                             get_id(k+1,i+1,j+1), get_id(k  ,i+1,j));
                mesh.add_tet(get_id(k  ,i+1,j+1), get_id(k  ,i+1,j),
                             get_id(k+1,i+1,j+1), get_id(k  ,i  ,j+1));
                mesh.add_tet(get_id(k+1,i  ,j+1), get_id(k+1,i+1,j+1),
                             get_id(k+1,i  ,j  ), get_id(k  ,i  ,j+1));
                mesh.add_tet(get_id(k  ,i  ,j  ), get_id(k  ,i  ,j+1),
                             get_id(k+1,i  ,j  ), get_id(k  ,i+1,j));
                mesh.add_tet(get_id(k+1,i+1,j+1), get_id(k  ,i+1,j),
                             get_id(k+1,i  ,j  ), get_id(k  ,i  ,j+1));
            }
        }

        return SUCC_RETURN;
    }

    static inline int get_id(int x, int y, int z)
    {
        return y*(NX+1)*(NZ+1) + z*(NX+1) + x;
    }
};
#endif
