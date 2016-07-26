/******************************************************************************
 *  File: FractureImpRecorder.h
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
#ifndef IO_FRACTURE_IMPULSE_RECORDER_H
#   define IO_FRACTURE_IMPULSE_RECORDER_H

#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include "StaticStressSolver.h"
#include "geometry/FixVtxTetMesh.hpp"
#include "rigid/LSCollisionRigidBody.hpp"


class FractureSimulator;

/*!
 * fracture 
 */
class FractureImpRecorder
{
    friend class FractureSimulator;

    public:
        typedef FixVtxTetMesh<REAL>                     TMesh;
        typedef LSCollisionRigidBody<REAL, TMesh>       TRigidBody;
        typedef FreeStaticStressSolver                  TStressSolver;

    private:
        struct ObjRec
        {
            TStressSolver*      psolver;
            std::vector<int>    surfIdMap;      // surface id map from surface vtx to tet vtx
            std::set<int>       appliedVtx;     // the vertex where impulses applied on

            ObjRec():psolver(NULL) {} 
            ObjRec(TStressSolver* ps):psolver(ps) 
            { 
                surfIdMap.clear();
                appliedVtx.clear();
            }
        };

    public:
        ~FractureImpRecorder()
        {   fout_.close(); }

        void init(const char* file, int precision = 12)
        {
            fout_.close();
            fout_.open(file);
            fout_ << std::setprecision(precision);
        }

        /*
         * record the impulse, applied on the vertex (id = cRec.vtxId)
         * of the rigid object (body).
         */
        void record_impulse(REAL ts, const Vector3<REAL>& imp, int vtxId, 
                            const TRigidBody* body, bool surfVtx);
        /*
         * record the impulse applied to the specified point
         */
        void record_impulse(REAL ts, const Vector3<REAL>& imp, 
                            const Point3<REAL>& pt, const TRigidBody* body);
        void record_extra_impulse(REAL ts, const Vector3<REAL>& imp, 
                                  int vtxId, const TRigidBody* body);

        const std::vector<int>& surface_id_map(int id) const
        {   
            std::map<int, ObjRec*>::const_iterator pt = idMap_.find(id);
            return pt->second->surfIdMap;
        }

        void add_rigid_body(int id, TStressSolver* psolver);

        /*
         * This method is called when a object is set to be unbreakable
         */
        void set_unbreakable(int id);

        void remove_rigid_body(int id)
        {   
            delete idMap_[id];
            idMap_.erase(id); 
        }

        /*
         * clean the impulses on each object
         */
        void time_step_begin();

        std::map<int, ObjRec*>  idMap_;      // id map from surface to tet mesh for each obj
    private:
        std::ofstream           fout_;
        REAL                    invStepSize_;
};

#endif
