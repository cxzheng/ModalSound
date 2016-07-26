/******************************************************************************
 *  File: DistFunc.hpp
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
#ifndef RIGID_DIST_FUNC_HPP
#   define RIGID_DIST_FUNC_HPP

/*
 * TriangleForCollision is the data structure for each triangle 
 * It maintains the info for transforming a triangle into the layout shown
 * below. 
 *
 *
 *    y^    \G E/
 *     |  N21\ /N12
 *     |      2
 *     | I   / \
 *     |    /   \  F
 * \N22|   /     \       /
 *  \  |  /S2   S1\     /N11
 *   \ | /    J    \   /   
 *  H \|/           \ /  D
 * ----0-------------1--------->
 *     |             |         x
 *     |             |
 *  A  |      C      |    B
 *     |             |
 */
template <typename T>
struct TriangleForCollision
{
    /*
     * vertex id of the triangle
     * The first vertex is gonna be transformed to origin
     */
    int         vtxId[3];

    /*
     * Transformation matrix, making v0, v1, v2 on X-Y plane
     * NOTE that, the matrix is constructed based on the INITIAL
     * state of v0, v1, and v2, which means the transformation
     * could be pre-computed.
     *
     * Q[0]: make v0-v1 along x-axis, v0 is at origin
     * Q[1]: make v1-v2 along x-axis, v1 is at origin
     * Q[2]: make v2-v0 along x-axis, v2 is at origin
     */
    Matrix3<T>  Q[3];

    /*
     * The inverse of Q
     */
    Matrix3<T>  QT[3];
    Vector3<T>  x0[3];

    T           sidelen[3];  // {|V1 - V0|, |V1-V2|, |V2-V0|
};

///////////////////////////////////////////////////////////////////////////////////////

/*
 * Helper class for computing the signed distance value from a point to a rigid body.
 */
template <typename T, class TCollProc>
class DistFunc
{
    private:
        TCollProc*      mp_cproc;

    public:
        DistFunc(TCollProc* pcproc):mp_cproc(pcproc) { }

        /* evaluate the level-set value at the given position 
         * The given position is in the rigid object's initial configuration
         */
        T operator() (const Point3<T>& pt)
        {
            _MinPtRec nearestPtRec;
            nearestPtRec.distSqr = 1E+100;

            if ( signed_distance(mp_cproc->bounding_tree_root(), pt, nearestPtRec) )
            {
                T ret = sqrt(nearestPtRec.distSqr);
                // determine the sign of the find shortest distance
                switch ( nearestPtRec.type )
                {
                    case 0:     // nearest point is on triangle
                        if ( (pt - nearestPtRec.closestPt).dotProduct(
                                    mp_cproc->m_tglPseudoNml[nearestPtRec.ids[0]]) < 0. )
                            return -ret;
                        else
                            return ret;
                    case 1:     // nearest point is on edge
                        if ( (pt - nearestPtRec.closestPt).dotProduct(
                                    mp_cproc->edge_pseudo_normal(
                                        nearestPtRec.ids[0], nearestPtRec.ids[1])) < 0. )
                            return -ret;
                        else
                            return ret;
                    case 2:     // nearest point is on vertex
                        if ( (pt - mp_cproc->mp_mesh->vertex(nearestPtRec.ids[0])).dotProduct(
                                    mp_cproc->m_vtxPseudoNml[nearestPtRec.ids[0]]) < 0. )
                            return -ret;
                        else
                            return ret;
                    default:
                        fprintf(stderr, "ERROR: Unknown _MinPtRec.type %d\n",  nearestPtRec.type);
                        exit(1);
                }
            }

            // should not reach here
            fprintf(stderr, "ERROR: fail to compute the distance field value\n");
            exit(1);
        }

    private:

        struct _MinPtRec
        {
            T   distSqr;
            int type;   // 0: nearest point is on triangle 
                        // 1: nearest point is on the edges
                        // 2: nearest point is on a vertex
            int ids[2]; // type = 0: ids[0] is the triangle id
                        // type = 1: ids[0],ids[1] is two ending point on that edge
                        // type = 2: ids[0] is the vertex id

            Point3<T>   closestPt;  // the closest point to the question position
                                    // if type = 3, this field is meaningless, because
                                    // the closest point is just the vertex specified
                                    // in ids[0]
        };

        /*
         * travel the OBB hierarchy to find the signed distance value
         * NOTE that pt should be in the rigid object's initial configuration
         *
         * return if the outRec is updated by this function call, if the given tree node
         *        is pruned out, false is returned
         */
        bool signed_distance(const typename TCollProc::TTreeNode* node, 
                const Point3<T>& pt, _MinPtRec& outRec)
        {
            // estimate the distance from the given pt to the bounding box of the current tree node
            if ( node->dist_sqr_to_node(pt) >= outRec.distSqr ) return false;

            if ( node->is_leaf() )
            {
                bool updated = false;
                // go over each triangle maintained on this tree node
                // compute the shortest distance to that triangle
                for(int cid = 0;cid < 2;++ cid)
                {
                    const std::vector<int>& ts = node->triangles(cid);
                    for(size_t tid = 0;tid < ts.size();++ tid)
                        updated |= shortest_distance(ts[tid], pt, outRec);
                }

                return updated;
            }

            bool bl = signed_distance(node->left_child(), pt, outRec);
            bool br = signed_distance(node->right_child(), pt, outRec);
            
            return (bl || br);
        }

        /*
         * compute the shortest distance from the given point to the 
         * given triangle tglId.
         *
         * NOTE that the field distSqr of the output outRec is specified when this method 
         *      is called. The input outRec.distSqr serves as an estimation about the 
         *      shortest distance so far. If the distance^2 to this triangle is larger 
         *      than that value, no change is made for outRec.
         *
         *    y^    \G E/
         *     |  N21\ /N12
         *     |      2
         *     | I   / \
         *     |    /   \  F
         * \N22|   /     \       /
         *  \  |  /S2   S1\     /N11
         *   \ | /    J    \   /   
         *  H \|/           \ /  D
         * ----0-------------1--------->
         *     |             |         x
         *     |             |
         *  A  |      C      |    B
         *     |             |
         *
         *     0        1    
         *    / \      / \   
         *   /   \    /   \  
         *  1-----2  2-----0 
         * 
         */
        bool shortest_distance(int tglId, const Point3<T>& pt, 
                _MinPtRec& outRec) const
        {
            const TriangleForCollision<T>& tglInfo = 
                mp_cproc->triangle_info_for_collision(tglId);

            T d0;
            Point3<T> P;

            for(int i = 0;i < 3;++ i)
            {
                // transform the point into this triangle's local coordinate system
                P = tglInfo.Q[i] * pt + tglInfo.x0[i];
                if ( !i ) d0 = P.z * P.z;

                if ( P.y > 0 ) continue;

                if ( P.x <= 0. )   // in region A
                {
                    const T tmpd = M_SQR(P.x) + M_SQR(P.y) + d0;
                    if ( tmpd < outRec.distSqr )
                    {
                        outRec.distSqr = tmpd;
                        outRec.type = 2;
                        outRec.ids[0] = tglInfo.vtxId[i];
                        return true;
                    }
                    return false;
                }
                else if ( P.x > tglInfo.sidelen[i] ) 
                {   // in region B
                    const T tmpd = M_SQR(P.y) + M_SQR(P.x - tglInfo.sidelen[i]) + d0;
                    if ( tmpd < outRec.distSqr )
                    {
                        outRec.distSqr = tmpd;
                        outRec.type = 2;
                        outRec.ids[0] = tglInfo.vtxId[(i+1)%3];
                        return true;
                    }
                    return false;
                }
                else
                {   // in region C
                    const T tmpd = M_SQR(P.y) + d0;
                    if ( tmpd < outRec.distSqr )
                    {
                        outRec.distSqr = tmpd;
                        outRec.type = 1;
                        outRec.ids[0] = tglInfo.vtxId[i];
                        outRec.ids[1] = tglInfo.vtxId[(i+1)%3];

                        // transform the closest point back in initial configuration 
                        // QT*(pn - x0)
                        outRec.closestPt.set(P.x, (T)0, (T)0);
                        outRec.closestPt -= tglInfo.x0[i];
                        outRec.closestPt = tglInfo.QT[i] * outRec.closestPt;

                        return true;
                    }
                    return false;
                }
            }

            // P is inside of the triangle
            if ( d0 < outRec.distSqr )
            {
                outRec.distSqr = d0;
                outRec.type = 0;
                outRec.ids[0] = tglId;

                // transform the closest point back in initial configuration 
                // P = tglInfo.Q[2] * pt + tglInfo.x0[2];
                // QT*(pn - x0)
                outRec.closestPt.set(P.x, P.y, (T)0);
                outRec.closestPt -= tglInfo.x0[2];
                outRec.closestPt  = tglInfo.QT[2] * outRec.closestPt;

                return true;
            }
            return false;
        }
};

#endif

