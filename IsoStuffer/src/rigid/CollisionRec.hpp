/******************************************************************************
 *  File: CollisionRec.hpp
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
#ifndef RIGID_COLLISION_REC_HPP
#   define RIGID_COLLISION_REC_HPP

#include "geometry/Point3.hpp"
#include "linearalgebra/Vector3.hpp"

/*!
 * CollisionRec is associated with each rigid body to record the information about collisions.
 * The information includes:
 * - the penetration depth at each vertex
 * - direction of collision impulse at the vertex
 */
template <typename T>
struct CollisionRec
{
    T           depth;      // negative value
    Point3<T>   pt;         // the collision point in predicted configuration
    Vector3<T>  impulseDir;
    Vector3<T>  vrel;
    T           vnrel;
    T           eps;        // restitution coefficient
    T           mu;         // friction coefficient
#ifdef USE_RECORDER
    int         vtxId;      // the vtx id in surface triangle mesh. the referred 
                            // vertex is the current deepest vertx that penetrates
                            // into the other objects
#endif
};

#endif
