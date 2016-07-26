/******************************************************************************
 *  File: RigidTetMesh.hpp
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
#ifndef GEOMETRY_RIGID_TET_MESH_HPP
#   define GEOMETRY_RIGID_TET_MESH_HPP

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T>
class RigidTetMesh
{
    public:
        typedef Tuple4ui    TetIdx;  // vertices' index for each tet vertex
};

#ifdef USE_NAMESPACE
}
#endif
#endif
