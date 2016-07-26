/******************************************************************************
 *  File: RigidObjDispRecorder.h
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
#ifndef IO_RIGID_DISP_RECORDER_H
#   define IO_RIGID_DISP_RECORDER_H

#include <fstream>
#include "config.h"
#include "geometry/FixVtxTetMesh.hpp"
#include "rigid/LSCollisionRigidBody.hpp"

/*
 * Record the displacement (translation and rotation into files)
 *
 * <d_ts> 
 * <i_id1>   # first rigid object 
 * <d_translate.x> <d_translate.y> <d_translate.z> 
 * <d_rot.w> <d_rot.x> <d_rot.y> <d_rot.z>
 * <i_id2>   # 2nd rigid obj
 * <d_translate.x> <d_translate.y> <d_translate.z> 
 * <d_rot.w> <d_rot.x> <d_rot.y> <d_rot.z>
 * <-1>      # label the ending of this timestep
 * ...
 * <d_ts>
 * <i_id1>   # first rigid object ...
 *
 */
class RigidObjDispRecorder
{
    public:
        typedef FixVtxTetMesh<REAL>                     TMesh;
        typedef LSCollisionRigidBody<REAL, TMesh>       TRigidBody;

        void init(const char* file)
        {
            m_fout.close();
            m_fout.open(file, std::ios::binary);
        }

        /*
         * Begin recording the displacement (both translation and rotation)
         * for the given time moment
         */
        inline void begin_record(double ts)
        {
            m_fout.write((const char*)&ts, sizeof(double));
        }

        /*
         * write the displacement of the given rigid body into file
         */
        void record_displacement(const TRigidBody* body)
        {
            const int id = body->id();
            m_fout.write((const char*)&id, sizeof(int));

            const Point3<REAL>& c = body->mass_center();
            m_fout.write((const char*)&c, sizeof(Point3<REAL>));
            
            const Quaternion<REAL>& r = body->rotation();
            m_fout.write((const char*)&r, sizeof(Quaternion<REAL>));
        }

        inline void end_record()
        {
            const int END = -1;
            m_fout.write((const char*)&END, sizeof(int));
            // flush data after each time step
            m_fout.flush();
        }

        ~RigidObjDispRecorder()
        {
            m_fout.close();
        }

    private:
        std::ofstream       m_fout;
};

#endif
