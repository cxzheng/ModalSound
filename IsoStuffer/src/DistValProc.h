/******************************************************************************
 *  File: DistValProc.h
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
#ifndef DIST_VALUE_PROC_H
#   define DIST_VALUE_PROC_H

#include "config.h"

#ifdef USE_HASH_MAP
#   include <unordered_map>
#else
#   include <tr1/unordered_map>
#endif

#include "geometry/TriangleMesh.hpp"
#include "geometry/OBBTree.hpp"
#include "generic/null_type.hpp"
#include "rigid/DistFunc.hpp"

class DistValProc
{
    template <typename T, class TCollProc> friend class DistFunc;

    public:
        typedef OBBTree< double, TriangleMesh<double>, carbine::NullType >  TOBBTree;
        typedef TOBBTree::TNode                                             TTreeNode;

        DistValProc(const TriangleMesh<double>* mesh, int res, int nlevel);


        const TTreeNode* bounding_tree_root() const
        {   return m_obbtree.root(); }

        const Vector3d& edge_pseudo_normal(int a, int b)
        {
            assert(a != b);

            if ( a > b ) std::swap(a, b);
            return m_edgePseudoNml[a][b];
        }
        
        const TriangleForCollision<double>& triangle_info_for_collision(int tid) const
        {  
            assert(tid < m_tglCollInfo.size() && tid >= 0);
            return m_tglCollInfo[tid];
        }

        const Tuple3i& resolution() const
        {   return m_res; }

        double grid_size() const
        {   return m_gds; }
        
        const Point3d& min_point() const
        {   return m_minCorner; }

        /*!
         * \param res    The lowest level resolution (the highest resolution) of the
         *               bounding box.
         * \param margin In implementation, the object is gonna be surrounded by a 
         *               bounding box where the tet is filled in. In order to detect
         *               the tet which is cross the object surface, the bounding box
         *               cannot be very tight to the object, there must be some margin.
         */
        void update_mesh_params(int res, int nlevel, int margin = 7);

    private:
        void compute_tgl_pseudo_normals();
        void compute_vtx_pseudo_normals();
        void compute_edge_pseudo_normals();
        void compute_tgl_collision_info();

    private:
        const TriangleMesh<double>*     mp_mesh;
        TOBBTree                        m_obbtree;
        double                          m_gds;          // grid size in the bounding volume
        Tuple3i                         m_res;
        Point3d                         m_minCorner;    // minimum corner of the bounding volume
        Point3d                         m_minP, m_maxP; // minimum/maximum point of the given triangle mesh
        
        // the pre-computed pseudo normal information, and triangle transformation info for computing
        // the signed distance field. 
        // (Implement the paper: Generating Signed Distance Fields From Triangle Meshes)
        // NOTE that all the information are computed at the rigid obj's initial configuration.
        std::vector< Vector3d >             m_tglPseudoNml; // pseudo normal at each triangle
        std::vector< Vector3d >             m_vtxPseudoNml; // pseudo normal on each vertex
#ifdef USE_UNORDERED_MAP
        std::tr1::unordered_map< int, std::tr1::unordered_map< int, Vector3d > >    m_edgePseudoNml;
#else
        std::unordered_map< int, std::unordered_map< int, Vector3d > >    m_edgePseudoNml;
#endif
        std::vector< TriangleForCollision<double> >     m_tglCollInfo;  // collision information for each triangle
};

#endif
