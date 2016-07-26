/******************************************************************************
 *  File: CollisionConstraint.hpp
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
#ifndef RIGID_COLLISION_CONSTRAINT_HPP
#   define RIGID_COLLISION_CONSTRAINT_HPP

#include "geometry/Point3.hpp"
#include "rigid/CollisionRec.hpp"

template <typename T, class TCollProc>
class CollisionConstraint
{
    public:
        typedef typename TCollProc::TTreeNode   TTreeNode;

    public:
        /* return if this constraint is fixed constraint 
         * the ground, static wall are all fixed constraints
         */
        virtual bool is_fixed() const { return false; }

        virtual bool deepest_collision(
                TCollProc*, TTreeNode*, CollisionRec<T>&)
        { return false; }

        virtual bool is_colliding(TCollProc*, TTreeNode*)
        { return false; }
};

/*
 * collision constraint for ground
 */
template <typename T, class TCollProc>
class GroundCollisionConstraint : public CollisionConstraint<T, TCollProc>
{
    public:
        typedef typename TCollProc::TTreeNode       TTreeNode;
        typedef typename TCollProc::TTreeNodeData   TTreeNodeData;

    public:
        GroundCollisionConstraint(T h = 0, T mu = 0.1, T eps = 0.9):
            m_gdHeight(h), m_mu(mu), m_eps(eps)
        { }

        bool is_fixed() const { return true; }

        bool deepest_collision(TCollProc* cp, TTreeNode* node, CollisionRec<T>& outRec)
        {
            // check if the current bounding box is colliding
            TTreeNodeData* data = node->data();
            if ( data->ts < cp->pred_timestamp() )
                cp->update_tree_node_state(data, node);
            const Tuple3<T> tr = node->r(); 

            T cy;
            int dx, dy, dz;
            T tx, ty, tz;
            for(dx = 0, tx = tr.x;dx < 2;++ dx, tx *= -1)
            for(dy = 0, ty = tr.y;dy < 2;++ dy, ty *= -1)
            for(dz = 0, tz = tr.z;dz < 2;++ dz, tz *= -1)
            {
                cy = data->predc0.y + data->predR.cols[0].y*tx +
                        data->predR.cols[1].y*ty + 
                        data->predR.cols[2].y*tz;
                if ( cy < m_gdHeight ) goto L_COLLIDING;
            }
            return false;

L_COLLIDING:
            if ( node->is_leaf() )
                return deepest_collision_vtx(cp, node, outRec);
            else
            {
                bool b1 = deepest_collision(cp, node->left_child(), outRec);
                bool b2 = deepest_collision(cp, node->right_child(), outRec);
                return b1 || b2;
            }
        }

        bool is_colliding(TCollProc* cp, TTreeNode* node)
        {
            // check if the current bounding box is colliding
            TTreeNodeData* data = node->data();
            if ( data->ts < cp->pred_timestamp() )
                cp->update_tree_node_state(data, node);
            const Tuple3<T> tr = node->r(); 

            T cy;
            int dx, dy, dz;
            T tx, ty, tz;
            for(dx = 0, tx = tr.x;dx < 2;++ dx, tx *= -1)
            for(dy = 0, ty = tr.y;dy < 2;++ dy, ty *= -1)
            for(dz = 0, tz = tr.z;dz < 2;++ dz, tz *= -1)
            {
                cy = data->predc0.y + data->predR.cols[0].y*tx +
                        data->predR.cols[1].y*ty + 
                        data->predR.cols[2].y*tz;
                if ( cy < m_gdHeight ) goto L_CHECK_COLLIDING;
            }
            return false;

L_CHECK_COLLIDING:
            if ( node->is_leaf() )
                return has_collision_vtx(cp, node);
            else
                return ( is_colliding(cp, node->left_child()) || 
                         is_colliding(cp, node->right_child()) );
        }

    private:
        bool has_collision_vtx(TCollProc* cp, TTreeNode* node) const
        {
            std::set<int> vids;
            const std::vector<Tuple3ui>& tvIds = cp->surface_mesh()->surface_indices();

            for(int cid = 0;cid < 2;++ cid)
            {
                const std::vector<int>& ts = node->triangles(cid);
                for(size_t i = 0;i < ts.size();++ i)    // iterate on triangle
                for(int vid = 0;vid < 3;++ vid)         // iterate on vertex
                {
                    if ( vids.count(tvIds[ts[i]][vid]) ) continue;
                    vids.insert(tvIds[ts[i]][vid]);

                    const Point3<T>& ppt = cp->predicted_surf_vtx_position(
                            tvIds[ts[i]][vid]);
                    if ( ppt.y < m_gdHeight ) return true;
                }
            }
            return false;
        }

        bool deepest_collision_vtx(TCollProc* cp, 
                TTreeNode* node, CollisionRec<T>& outRec) const
        {
            bool ret = false;
            std::set<int> vids;
            const std::vector<Tuple3ui>& tvIds = cp->surface_mesh()->surface_indices();

            for(int cid = 0;cid < 2;++ cid)
            {
                const std::vector<int>& ts = node->triangles(cid);
                for(size_t i = 0;i < ts.size();++ i)
                for(int vid = 0;vid < 3;++ vid)
                {
                    if ( vids.count(tvIds[ts[i]][vid]) ) continue;
                    vids.insert(tvIds[ts[i]][vid]);

                    const Point3<T>& ppt = cp->predicted_surf_vtx_position(
                            tvIds[ts[i]][vid]);
                    if ( ppt.y < m_gdHeight && outRec.depth > ppt.y - m_gdHeight)
                    {
                        const Vector3<T> vnow = cp->rigid_body()->predicted_velocity(ppt); 
                        //const Vector3<T> vnow = cp->rigid_body()->current_velocity(ppt); 
                        if ( vnow.y < 0 )
                        {
                            outRec.depth = ppt.y - m_gdHeight;
                            outRec.pt    = ppt;
                            outRec.impulseDir.set(0., 1., 0.);
                            outRec.vnrel = vnow.y;
                            outRec.vrel  = vnow;
                            outRec.eps   = fmin(m_eps, cp->rigid_body()->rest_coeff());
                            outRec.mu    = fmax(m_mu, cp->rigid_body()->friction_coeff());
#ifdef USE_RECORDER
                            outRec.vtxId = tvIds[ts[i]][vid];
#endif
                            ret = true;
                        }
                    }
                }
            }
            return ret;
        }

    private:
        T           m_gdHeight;

        T           m_mu;       // friction coefficient
        T           m_eps;      // restitution coefficient
};

#endif
