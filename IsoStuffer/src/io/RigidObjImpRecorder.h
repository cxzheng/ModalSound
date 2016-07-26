/******************************************************************************
 *  File: RigidObjImpRecorder.h
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
#ifndef IO_RIGID_IMPULSE_RECORDER_H
#   define IO_RIGID_IMPULSE_RECORDER_H

#include <fstream>
#include "geometry/FixVtxTetMesh.hpp"
#include "rigid/LSCollisionRigidBody.hpp"

/*
 * Record impulses applied on objects into file. These impulses 
 * will be used to drive the vibration later when generating the 
 * sound.
 *
 * The recording in the file is gonna to be
 * <time>    <obj id>    <applied vtx id>    <impulse>
 * ...
 *
 * NOTE that the "applied vtx id" is the ID on surface triangle mesh
 *
 */
class RigidObjImpRecorder
{
    public:
        typedef FixVtxTetMesh<REAL>                     TMesh;
        typedef LSCollisionRigidBody<REAL, TMesh>       TRigidBody;

    public:
        ~RigidObjImpRecorder()
        {
            m_fout.close();
        }

        void init(const char* file, int precision = 12);
        /*
         * record the impulse, applied on the vertex (id = cRec.vtxId)
         * of the rigid object (body).
         */
        void record_impulse(REAL ts, const Vector3<REAL>& imp, 
                int vtxId, const TRigidBody* body, bool surfVtx);
        void record_impulse(REAL ts, const Vector3<REAL>& imp, 
                const Point3<REAL>& pt, const TRigidBody* body);

    private:
        std::ofstream       m_fout;
};

#endif
