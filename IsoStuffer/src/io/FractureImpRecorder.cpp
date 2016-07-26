/******************************************************************************
 *  File: FractureImpRecorder.cpp
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
#include "FractureImpRecorder.h"
#include "utils/print_msg.h"

using namespace std;

static const REAL IMPULSE_SCALE = 1800.; //1./0.0003;

void FractureImpRecorder::add_rigid_body(int id, TStressSolver* psolver)
{   
    assert(!idMap_.count(id));
    ObjRec * ptr;
    idMap_[id] = ptr = new ObjRec(psolver);

    if ( !psolver ) return;

    printf("INFO: compute surface ID map\n");
    map<int, int> idmap;    // map from tet id to surface id
    psolver->mesh()->surface_id_map(idmap);

    // map from surface id to tet id
    ptr->surfIdMap.resize(idmap.size());
    const map<int, int>::iterator end = idmap.end();
    for(map<int, int>::iterator it = idmap.begin();it != end;++ it)
        ptr->surfIdMap[it->second] = it->first;
}

void FractureImpRecorder::set_unbreakable(int id)
{
    if ( !idMap_.count(id) )
    {
        PRINT_WARNING("No obj id=%d in the FractureImpRecorder table\n");
        return;
    }
    ObjRec* ptr = idMap_[id];
    ptr->psolver = NULL;
    ptr->surfIdMap.clear();
    ptr->appliedVtx.clear();
}

/*
 * NOTE: Make sure here the input vector imp is already in object's 
 *       rest(initial) configuration
 *       Also make sure the input vtxId is the vertex ID in tet 
 *       mesh
 */
void FractureImpRecorder::record_extra_impulse(REAL ts, 
        const Vector3<REAL>& imp, int vtxId, const TRigidBody* body)
{
    fout_ << ts << ' ' << body->id() << ' ' 
          << vtxId << ' '      // vtxId is 0-based
          << imp.x << ' '
          << imp.y << ' '
          << imp.z << " T" << std::endl; 
}

void FractureImpRecorder::record_impulse(REAL ts, const Vector3<REAL>& imp,
        int vtxId, const TRigidBody* body, bool surfVtx)
{
    // map the impulse vector to object's rest configuration
    const Vector3<REAL> impVec = body->predicted_inverse_rotation().rotate(imp);

    ObjRec* objrec = idMap_[body->id()];

    // for unbreakable objects
    if ( !objrec->psolver )
    {
        fout_ << ts << ' ' << body->id() << ' ' 
              << vtxId << ' '      // vtxId is 0-based
              << impVec.x << ' '
              << impVec.y << ' '
              << impVec.z << ' '
              << (surfVtx ? 'S' : 'T') << std::endl; 
        return;
    }

    if ( surfVtx ) vtxId = objrec->surfIdMap[vtxId];
    int numFixed = body->mesh()->num_fixed_vertices();
    if ( vtxId < numFixed ) return; // ignore the impulses applied on fixed vertices

    objrec->appliedVtx.insert(vtxId);
    vector< Vector3<REAL> >& fs = objrec->psolver->current_impulses();
#ifdef WITH_QUASI_STATIC_OBJS
    fs[vtxId].scaleAdd(IMPULSE_SCALE, impVec);
#else
    fs[vtxId - numFixed].scaleAdd(IMPULSE_SCALE, impVec);
#endif

    //cerr << "IMPULSES: " << imp << impVec << body->predicted_inverse_rotation() 
    //     << body->predicted_rotation() << endl;

    fout_ << ts << ' ' << body->id() << ' ' 
          << vtxId << ' '      // vtxId is 0-based
          << impVec.x << ' '
          << impVec.y << ' '
          << impVec.z << " T" << std::endl; 
}

void FractureImpRecorder::record_impulse(REAL ts, const Vector3<REAL>& imp,
        const Point3<REAL>& pt, const TRigidBody* body)
{
#ifdef USE_RECORDER
    // map the impulse vector to object's rest configuration
    const Vector3<REAL> impVec = body->predicted_inverse_rotation().rotate(imp);

    //// transform \pt from "predicted" configuration into rest configuration
    const Point3<REAL> ipt = body->initial_predicted_position(pt);
    //// find the nearest vertex from \ipt
    //   vid is the id in surface mesh because the kd-tree is built only on surface mesh
    int vtxId = body->kdtree().find_nearest(ipt);
    ObjRec* objrec = idMap_[body->id()];

    if ( !objrec->psolver )
    {
        fout_ << ts << ' ' << body->id() << ' ' 
              << vtxId << ' '      // vtxId is 0-based
              << impVec.x << ' '
              << impVec.y << ' '
              << impVec.z << " S" << std::endl;
        return;
    }

    vtxId = objrec->surfIdMap[vtxId];
    int numFixed = body->mesh()->num_fixed_vertices();
    if ( vtxId < numFixed ) return; // ignore the impulses applied on fixed vertices

    objrec->appliedVtx.insert(vtxId);
    vector< Vector3<REAL> >& fs = objrec->psolver->current_impulses();
#ifdef WITH_QUASI_STATIC_OBJS
    fs[vtxId].scaleAdd(IMPULSE_SCALE, impVec);
#else
    fs[vtxId - numFixed].scaleAdd(IMPULSE_SCALE, impVec);
#endif

    fout_ << ts << ' ' << body->id() << ' '
          << vtxId << ' '
          << impVec.x << ' '
          << impVec.y << ' '
          << impVec.z << " T" << std::endl;
#endif
}

/*!
 * clean the current impulses
 */
void FractureImpRecorder::time_step_begin()
{
    const map<int, ObjRec*>::iterator end = idMap_.end();
    for(map<int, ObjRec*>::iterator it = idMap_.begin();it != end;++ it)
    {
        if ( !it->second->psolver ) continue;

        it->second->appliedVtx.clear();
        vector< Vector3<REAL> >& imp = it->second->psolver->current_impulses();
        memset(&imp[0], 0, sizeof(Vector3<REAL>)*imp.size());
    }
}

