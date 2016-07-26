/******************************************************************************
 *  File: SurfForcedTetMesh.hpp
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
#ifndef SURFACE_FORCED_TET_MESH_HPP
#   define SURFACE_FORCED_TET_MESH_HPP

#include <assert.h>
#include <vector>
#include "linearalgebra/Vector3.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T, class TMesh>
struct SurfForcedTetMesh
{
    /* the referred tet mesh */
    TMesh*      pmesh;

    /*
     * map each tet to the tet id in the original un-fractured object
     */
    std::vector<int>    origTids;

    SurfForcedTetMesh(TMesh* mesh):pmesh(mesh)
    {
        assert(pmesh != NULL);
    }

    SurfForcedTetMesh():pmesh(NULL) 
    { }

    void init_tet_map()
    {
        origTids.resize(pmesh->num_tets());
        for(size_t i = 0;i < origTids.size();++ i)
            origTids[i] = i;
    }
};

    /* 
     * the triangles indexed by its three vertices where there are stress
     * force applied on
    TForcedTgl   forcedTgls;
     */
#ifdef USE_NAMESPACE
}
#endif
#endif
