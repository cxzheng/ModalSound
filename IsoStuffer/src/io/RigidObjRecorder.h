/******************************************************************************
 *  File: RigidObjRecorder.h
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
#ifndef IO_RIGID_OBJ_RECORDER_H
#   define IO_RIGID_OBJ_RECORDER_H

#include "geometry/FixVtxTetMesh.hpp"
#include "deformable/stvk.h"
#include "rigid/LSCollisionRigidBody.hpp"
#include "linearalgebra/PardisoMatrix.hpp"

/*
 * save_geometry, save_mass_mat, save_stiff_mass are used to save the
 * necessary information for modal analysis.
 */
namespace RigidObjRecorder
{
    typedef FixVtxTetMesh<REAL>                     TMesh;
    typedef StVKMaterial                            TMaterial;
    typedef LSCollisionRigidBody<REAL, TMesh>       TRigidBody;

    /*!
     * save the surface mesh into an OBJ mesh file
     */
    int save_surface_mesh(const TRigidBody* body, const char* file);
    
    int save_geometry(const TRigidBody* body, const char* file);
    int save_mass_mat(const TRigidBody* body, const char* file);
    int save_stiff_mat(const TRigidBody* body, const TMaterial* pmaterial, 
                       const char* file);
    int save_displacement(double ts, const TetMesh<double>* mesh, 
                          const char* file);
}

#endif
