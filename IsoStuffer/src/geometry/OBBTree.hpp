/******************************************************************************
 *  File: OBBTree.hpp
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
#ifndef OBB_TREE_HPP
#   define OBB_TREE_HPP

#include <vector>
#include "OBBTreeNode.hpp"
#include "geometry/Triangle.hpp"
#include "utils/math.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

/*!
 * Implement the oriented bounding box
 */
template <typename T, class _TMesh, class _TData>
class OBBTree
{
    public:
        typedef OBBTree<T, _TMesh, _TData>      TSelf;
        typedef OBBTreeNode<T, TSelf>           TNode;
        typedef _TMesh                          TMesh;
        typedef _TData                          TData;

        struct moment
        {
            T           area;
            Point3<T>   centroid;
            Matrix3<T>  C;
        };

        static const int    LEAF_SIZE;  // threshold for the number of triangles maintained on a leaf

    public:
        OBBTree(const TMesh* pmesh)     // initalize with the given triangle mesh
        {  init(pmesh); }

        ~OBBTree();

        void init(const TMesh* pmesh);

        const TNode* root() const
        {  return mp_root; }

        TNode* root()
        {  return mp_root; }

        const std::vector<moment>& triangle_moments() const
        {  return m_moments; }

    private:
        void update_moments();

    private:
        TNode*                      mp_root;
        std::vector<moment>         m_moments;
};

///////////////////////////////////////////////////////////////////////////////
template <typename T, class _TMesh, class _TData>
const int OBBTree<T, _TMesh, _TData>::LEAF_SIZE = 4;

template <typename T, class _TMesh, class _TData>
void OBBTree<T, _TMesh, _TData>::init(const TMesh* pmesh)
{
    if ( pmesh == NULL ) 
    {
        mp_root = NULL;
        return;
    }

    const std::vector< Tuple3ui >&  indices = pmesh->surface_indices();
    const std::vector< Point3<T> >& vtx     = pmesh->vertices();
    std::vector<int> tgls(indices.size());
    m_moments.resize(indices.size());

    const REAL TS = 1./3.;
    const REAL TS2 = 1./12.;

    //// compute moment for each triangle
    for(size_t tid = 0;tid < tgls.size();++ tid)
    {
        tgls[tid] = tid;

        const Point3<T>& v0 = vtx[indices[tid][0]]; 
        const Point3<T>& v1 = vtx[indices[tid][1]]; 
        const Point3<T>& v2 = vtx[indices[tid][2]];
        m_moments[tid].area = Triangle<T>::area(v0, v1, v2);
        m_moments[tid].centroid = (v0 + v1 + v2) * TS;

        if ( m_moments[tid].area < 1E-18 )
        {
            for(int k = 0;k < 3;++ k)
            {
                m_moments[tid].C.cols[k][k] = 
                    M_SQR(v0[k]) + M_SQR(v1[k]) + M_SQR(v2[k]);
                for(int j = k+1;j < 3;++ j)
                    m_moments[tid].C.cols[k][j] = m_moments[tid].C.cols[j][k] = 
                        v0[k]*v0[j] + v1[k]*v1[j] + v2[k]*v2[j];
            }
        }
        else
        {
            for(int k = 0;k < 3;++ k)
            {
                m_moments[tid].C.cols[k][k] = m_moments[tid].area * TS2 * 
                    (9.*M_SQR(m_moments[tid].centroid[k]) + 
                     M_SQR(v0[k]) + M_SQR(v1[k]) + M_SQR(v2[k]));
                for(int j = k+1;j < 3;++ j)
                {
                    m_moments[tid].C.cols[k][j] = m_moments[tid].C.cols[j][k] = 
                        m_moments[tid].area * TS2 * 
                        (9.*m_moments[tid].centroid[k]*m_moments[tid].centroid[j] +
                         v0[k]*v0[j] + v1[k]*v1[j] + v2[k]*v2[j]);
                }
            }
        }
    }  // end for 

    mp_root = new OBBTreeNode<T, TSelf>(this, pmesh, tgls, 0);
}

template <typename T, class _TMesh, class _TData>
OBBTree<T, _TMesh, _TData>::~OBBTree()
{
    delete mp_root;
}

#ifdef USE_NAMESPACE
}
#endif
#endif
