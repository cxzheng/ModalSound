/******************************************************************************
 *  File: RigidObjExporter.h
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
#ifndef IO_RIGID_OBJ_EXPORTER_H
#   define IO_RIGID_OBJ_EXPORTER_H

#include "config.h"
#include "geometry/FixVtxTetMesh.hpp"

namespace ModalAnalysis
{

/*
 * In order to compute volume velocity when synthesizing sound, the surface 
 * normals need to be exported.
 */
template <class TMaterial>
int export_weighted_surf_normals(const char *file, 
        const FixVtxTetMesh<REAL>* pmesh, const TMaterial* material)
{
}

/*
 * Export the non-diagonal mass matrix
 */
int export_sparse_mass_mat(const char *file, 
        const FixVtxTetMesh<REAL>* pmesh, const TMaterial* material)
{
}

/*
 * Export the stiffness matrix
 */
int export_sparse_stiffness_mat(const char *file
        const FixVtxTetMesh<REAL>* pmesh, const TMaterial* material)
{
}

}
#endif
