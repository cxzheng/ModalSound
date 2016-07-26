/******************************************************************************
 *  File: OBBTreeNode.hpp
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
#ifndef OBB_TREE_NODE_HPP
#   define OBB_TREE_NODE_HPP

#include <vector>
#include <set>
#include <algorithm>
#include "geometry/Point3.hpp"
#include "linearalgebra/Matrix3.hpp"
#include "linearalgebra/eig3.h"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T, class TTree>
class OBBTreeNode
{
    public:
        typedef typename TTree::TMesh   TMesh;
        typedef typename TTree::TData   TData;
        typedef OBBTreeNode<T, TTree>   TSelf;

    public:
        OBBTreeNode(TTree* ptree, 
                const TMesh* pmesh, 
                const std::vector<int>& tgls, 
                int level);
        ~OBBTreeNode();

        bool is_leaf() const
        {  return m_isLeaf; }
        int level() const
        {  return m_level; }
        inline int num_children() const
        {  return 2; }
        const TSelf** children() const
        {  return (const TSelf **)m_child; }
        TSelf** children() 
        {  return m_child; }
        const TSelf* left_child() const
        {  return m_child[0]; }
        TSelf* left_child() 
        {  return m_child[0]; }
        const TSelf* right_child() const
        {  return m_child[1]; }
        TSelf* right_child() 
        {  return m_child[1]; }

        const Point3<T>& c() const      // center of bounding box
        {  return m_c; }
        const Tuple3<T>& r() const      // radius in 3 dimensions
        {  return m_d; }
        const Matrix3<T>& R() const     // vector in 3 principle directions
        {  return m_R; }
        size_t size() const
        {  return m_tgls[0].size() + m_tgls[1].size(); }

        const std::vector<int>& triangles(int c) const
        {
            assert(c==0 || c==1);
            return m_tgls[c];
        }
        // NOTE: only available for leaf node
        const std::set<int>& vertices() const
        {   return m_vtx; }

        /*
         * like the same method in KD-tree, this method compute the distance
         * from the given point to the bouning box of this node.
         *
         * The returned value is also the lower bound of the distances from
         * the given point to the vertices inside this tree node.
         */
        T dist_sqr_to_node(const Point3<T>& pt) const;

        TData* data() {  return &m_data; }
        const TData* data() const {  return &m_data; }

    private:
        void eigen_decomp(const Matrix3<T>& CM);
        /*! build the bounding box */
        void build_bounding_box(const TMesh* pmesh,
                const Point3<T>& meanPt, 
                const std::vector<typename TTree::moment>& moments,
                const std::vector<int>& tgls);

    private:
        int         m_level;// level on the tree
        Point3<T>   m_c;    // center of the bounding box
        Matrix3<T>  m_R;    // the three axises for this bounding box in world-space coordinate
        Tuple3<T>   m_d;    // radius; half measure of a side length

        bool                m_isLeaf;
        TSelf*              m_child[2];
        std::vector<int>    m_tgls[2];  // the contained triangles, given by indices
        std::set<int>       m_vtx;      // the contained vertices, given by indices
                                        // only non-empty when this node is a leaf
        TData               m_data;
};

///////////////////////////////////////////////////////////////////////////////

/*!
 * - Compute the covariance matrix
 * - eigenvalue decomposition
 */
template <typename T, class TTree>
OBBTreeNode<T, TTree>::OBBTreeNode(TTree* ptree, const TMesh* pmesh, 
        const std::vector<int>& tgls, int level):m_level(level)
{
    //// compute the covariance matrix
    const std::vector<typename TTree::moment>& moments = ptree->triangle_moments();

    T           atot = 0;
    Point3<T>   wPtSum;
    Matrix3<T>  CM;
    for(size_t i = 0;i < tgls.size();++ i)
    {
        atot += moments[tgls[i]].area;
        wPtSum.scaleAdd(moments[tgls[i]].area, moments[tgls[i]].centroid);
        CM   += moments[tgls[i]].C;
    }

    const T invA = (T)1 / atot;
    for(int i = 0;i < 3;++ i)
    for(int j = 0;j < 3;++ j)
        CM.cols[i][j] -= wPtSum[j]*wPtSum[i]*invA;

    //// do eigen-decomposition
    // compute the eigenvectors, and sort them according to the eigenvalues
    eigen_decomp(CM);   // set the m_R matrix
    wPtSum *= invA;     // the mean point

    build_bounding_box(pmesh, wPtSum, moments, tgls);
    if ( (int)tgls.size() < TTree::LEAF_SIZE || 
            m_tgls[0].empty() || m_tgls[1].empty() )
    {
        //// create a leaf node
        m_child[0] = m_child[1] = NULL;
        m_isLeaf = true;
        // store all the vertices
        const std::vector<Tuple3ui>& indices = pmesh->surface_indices();
        for(size_t tid = 0;tid < tgls.size();++ tid)
        for(int i = 0;i < 3;++ i)
            m_vtx.insert(indices[tgls[tid]][i]);
    }
    else
    {
        //// create left child and right child
        m_isLeaf = false;
        m_child[0] = new TSelf(ptree, pmesh, m_tgls[0], level+1);
        m_child[1] = new TSelf(ptree, pmesh, m_tgls[1], level+1);
    }
}

template <typename T, class TTree>
OBBTreeNode<T, TTree>::~OBBTreeNode()
{
    delete m_child[0];
    delete m_child[1];
}

/*
 * Build the bounding box
 * - set m_d radius in 3 dimensions
 * - set m_c center of the bounding box
 */
template <typename T, class TTree>
void OBBTreeNode<T, TTree>::build_bounding_box(const TMesh* pmesh, 
        const Point3<T>& meanPt, 
        const std::vector<typename TTree::moment>& moments,
        const std::vector<int>& tgls)
{
    m_tgls[0].clear();
    m_tgls[1].clear();

    const std::vector<Tuple3ui>&      indices = pmesh->surface_indices();
    const std::vector< Point3<T> >&   vtx     = pmesh->vertices();
    REAL axdmp = m_R.cols[0].dotProduct(meanPt);
    Tuple3<REAL> minval(1E+100, 1E+100, 1E+100), 
                 maxval(-1E+100, -1E+100, -1E+100);
    /*
     * - min/max value along three principle direction 
     */
    for(size_t tid = 0;tid < tgls.size();++ tid)
    {
        // TODO: use cblas for better performance
        for(int vid = 0;vid < 3;++ vid)  // go through each vertex on this triangle
        for(int i = 0;i < 3;++ i)        // each dimension
        {
            REAL td = m_R.cols[i].dotProduct(vtx[indices[tgls[tid]][vid]]);
            minval[i] = fmin(minval[i], td);
            maxval[i] = fmax(maxval[i], td);
        }

        // use the center of the triangle to determine which child node it belongs to
        REAL td = m_R.cols[0].dotProduct(moments[tgls[tid]].centroid);
        if ( td <= axdmp )  // add to 1st child
            m_tgls[0].push_back(tgls[tid]);
        else
            m_tgls[1].push_back(tgls[tid]);
    }

    const Vector3<REAL> cc = (minval + maxval)*0.5;
    m_c = m_R * cc;
    m_d = (maxval - minval) * 0.5;
    m_d += 1E-10;
}

/*
 * - eigen decomposition
 * - sort it based on eigenvalues
 */
template <typename T, class TTree>
void OBBTreeNode<T, TTree>::eigen_decomp(const Matrix3<T>& CM)
{
    REAL V[3][3], d[3];
    eigen_decomposition((const double(*)[3])&CM, V, d);

    int maxId = 0;
    if ( fabs(d[1]) > fabs(d[maxId]) ) maxId = 1;
    if ( fabs(d[2]) > fabs(d[maxId]) ) maxId = 2;

    if ( maxId != 0 )
    {
        //// swap the maximum eigenvalue to be the first one
        std::swap(V[0][0], V[0][maxId]);
        std::swap(V[1][0], V[1][maxId]);
        std::swap(V[2][0], V[2][maxId]);
    }
    m_R.set(V[0][0], V[0][1], V[0][2],
            V[1][0], V[1][1], V[1][2],
            V[2][0], V[2][1], V[2][2]);
    assert((m_R.cols[0].normSqr()-1.) < 1E-9 &&
           (m_R.cols[1].normSqr()-1.) < 1E-9 &&
           (m_R.cols[2].normSqr()-1.) < 1E-9);
}

template <typename T, class TTree>
T OBBTreeNode<T, TTree>::dist_sqr_to_node(const Point3<T>& pt) const
{
    // coordinate in the local coordinate system of this node
    const Vector3<T> v = pt - m_c;
    T ret = 0;

    Vector3<T> coord(m_R.cols[0].dotProduct(v),
                     m_R.cols[1].dotProduct(v),
                     m_R.cols[2].dotProduct(v));
    for(int i = 0;i < 3;++ i)
    {
        if ( coord[i] < -m_d[i] ) 
            ret += M_SQR(coord[i] + m_d[i]);
        else if ( coord[i] >  m_d[i] )
            ret += M_SQR(coord[i] - m_d[i]);
    }
    return ret;
}

#ifdef USE_NAMESPACE
}
#endif
#endif
