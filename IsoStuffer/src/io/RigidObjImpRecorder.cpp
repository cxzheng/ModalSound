/******************************************************************************
 *  File: RigidObjImpRecorder.cpp
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
#include "RigidObjImpRecorder.h"
#include <iomanip>

using namespace std;

void RigidObjImpRecorder::init(const char* file, int precision)
{
    m_fout.close();
    m_fout.open(file);
    m_fout << setprecision(precision);
}

/*
 * - Transform the impulse vector from the current prediction to
 *   the object's initial(rest) configuration
 *
 * the last field in each line is a char, either T or S, indicating
 * if the given vertex id is the surface vtx id or tet mesh vtx id.
 * 'S' means surface vtx id, 'T' means tet id
 * The recording in the file is gonna to be
 * <time>    <obj id>    <applied vtx id>    <impulse>  <T/S>
 * ...
 */
void RigidObjImpRecorder::record_impulse(
        REAL ts, const Vector3<REAL>& imp, int vtxId,
        const TRigidBody* body, bool surfVtx)
{
    // map the impulse vector to object's rest configuration
    const Vector3<REAL> impVec = body->predicted_inverse_rotation().rotate(imp);

    m_fout << ts << ' ' << body->id() << ' ' 
           << vtxId << ' '      // vtxId is 0-based
           << impVec.x << ' '
           << impVec.y << ' '
           << impVec.z << ' ' 
           << (surfVtx ? 'S' : 'T') << std::endl;
}

/*
 * \pt is the point where the impulse is applied in the current "predicted" 
 * configuration. This method looks up the nearest vertex to \pt in the
 * vertices of the given rigid body \body, and records the impulse at that 
 * vertex
 */
void RigidObjImpRecorder::record_impulse(REAL ts, 
        const Vector3<REAL>& imp, const Point3<REAL>& pt,
        const TRigidBody* body)
{
#ifdef USE_RECORDER
    const Vector3<REAL> impVec = body->predicted_inverse_rotation().rotate(imp);
    //// transform \pt from "predicted" configuration into rest configuration
    const Point3<REAL> ipt = body->initial_predicted_position(pt);
    //// find the nearest vertex from \ipt
    //   vid is the id in tet mesh because the kd-tree is built from the entire 
    //   tet vertices
    int vid = body->kdtree().find_nearest(ipt);

    m_fout << ts << ' ' << body->id() << ' '
           << vid << ' '
           << impVec.x << ' '
           << impVec.y << ' '
           << impVec.z << " T" << std::endl;
#endif
}

