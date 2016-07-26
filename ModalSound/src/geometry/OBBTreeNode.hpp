#ifndef OBB_TREE_NODE_HPP
#   define OBB_TREE_NODE_HPP

#include <vector>
#include <algorithm>
#include <limits>
#include "generic/precision_type.hpp"
#include "Point3.hpp"
#include "sc/Matrix3.hpp"
#include "sc/eig3.h"
#include "utils/math.hpp"

#ifdef USE_NAMESPACE
namespace sploosh
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
        {  return isLeaf_; }
        int level() const
        {  return level_; }
        const TSelf** children() const
        {  return (const TSelf **)child_; }
        TSelf** children() 
        {  return child_; }
        const TSelf* left_child() const
        {  return child_[0]; }
        TSelf* left_child() 
        {  return child_[0]; }
        const TSelf* right_child() const
        {  return child_[1]; }
        TSelf* right_child() 
        {  return child_[1]; }

        const Point3<T>& c() const      // center of bounding box
        {  return c_; }
        const Tuple3<T>& r() const      // radius in 3 dimensions
        {  return d_; }
        const Matrix3<T>& R() const     // vector in 3 principle directions
        {  return R_; }
        size_t size() const
        {  return tgls_[0].size() + tgls_[1].size(); }

        const std::vector<int>& triangles(int c) const
        {
            assert(c==0 || c==1);
            return tgls_[c];
        }

        int is_disjoint(const TSelf* nb) const;

        /*
         * like the same method in KD-tree, this method compute the distance
         * from the given point to the bouning box of this node.
         *
         * The returned value is also the lower bound of the distances from
         * the given point to the vertices inside this tree node.
         */
        T dist_sqr_to_node(const Point3<T>& pt) const;

        TData* data() {  return &data_; }
        const TData* data() const {  return &data_; }

    private:
        void eigen_decomp(const Matrix3<T>& CM);
        /*! Build the bounding box */
        void build_bounding_box(const TMesh* pmesh,
                const Point3<T>& meanPt, 
                const std::vector<typename TTree::moment>& moments,
                const std::vector<int>& tgls);

    private:
        int         level_; // level on the tree
        Point3<T>   c_;     // center of the bounding box
        Matrix3<T>  R_;     // the three axises for this bounding box in world-space coordinate
        Tuple3<T>   d_;     // radius; half measure of a side length

        bool                isLeaf_;
        TSelf*              child_[2];
        std::vector<int>    tgls_[2];  // the contained triangles, given by indices
        TData               data_;
};

///////////////////////////////////////////////////////////////////////////////

/*!
 * - Compute the covariance matrix
 * - eigenvalue decomposition
 */
template <typename T, class TTree>
OBBTreeNode<T, TTree>::OBBTreeNode(TTree* ptree, const TMesh* pmesh, 
        const std::vector<int>& tgls, int level):level_(level)
{
    //// compute the covariance matrix
    const std::vector<typename TTree::moment>& moments = ptree->triangle_moments();

    T           atot = 0;
    Point3<T>   wPtSum;
    Matrix3<T>  CM;
    for(size_t i = 0;i < tgls.size();++ i)
    {
        atot += moments[tgls[i]].area;
        wPtSum.scale_add(moments[tgls[i]].area, moments[tgls[i]].centroid);
        CM   += moments[tgls[i]].C;
    }

    const T invA = (T)1 / atot;
    for(int i = 0;i < 3;++ i)
    for(int j = 0;j < 3;++ j)
        CM.cols[i][j] -= wPtSum[j]*wPtSum[i]*invA;

    //// do eigen-decomposition
    // compute the eigenvectors, and sort them according to the eigenvalues
    eigen_decomp(CM);   // set the R_ matrix
    wPtSum *= invA;     // the mean point

    build_bounding_box(pmesh, wPtSum, moments, tgls);
    if ( (int)tgls.size() < TTree::LEAF_SIZE || 
            tgls_[0].empty() || tgls_[1].empty() )
    {
        //// create a leaf node
        child_[0] = child_[1] = NULL;
        isLeaf_ = true;
    }
    else
    {
        //// create left child and right child
        isLeaf_ = false;
        child_[0] = new TSelf(ptree, pmesh, tgls_[0], level+1);
        child_[1] = new TSelf(ptree, pmesh, tgls_[1], level+1);
    }
}

template <typename T, class TTree>
OBBTreeNode<T, TTree>::~OBBTreeNode()
{
    delete child_[0];
    delete child_[1];
}

/*
 * Build the bounding box
 * - set d_ radius in 3 dimensions
 * - set c_ center of the bounding box
 */
template <typename T, class TTree>
void OBBTreeNode<T, TTree>::build_bounding_box(const TMesh* pmesh, 
        const Point3<T>& meanPt, 
        const std::vector<typename TTree::moment>& moments,
        const std::vector<int>& tgls)
{
    tgls_[0].clear();
    tgls_[1].clear();

    const std::vector<Tuple3ui>&      indices = pmesh->surface_indices();
    const std::vector< Point3<T> >&   vtx     = pmesh->vertices();

    const T inf = std::numeric_limits<T>::infinity();
    Tuple3<T> minval( inf,  inf,  inf),
              maxval(-inf, -inf, -inf);

    T axdmp = R_.cols[0].dot(meanPt);

    /*
     * - min/max value along three principle direction 
     */
    for(size_t tid = 0;tid < tgls.size();++ tid)
    {
        // TODO: use cblas for better performance
        for(int vid = 0;vid < 3;++ vid)  // go through each vertex on this triangle
        for(int i = 0;i < 3;++ i)        // each dimension
        {
            T td = R_.cols[i].dot(vtx[indices[tgls[tid]][vid]]);
            minval[i] = fmin(minval[i], td);
            maxval[i] = fmax(maxval[i], td);
        }

        // use the center of the triangle to determine which child node it belongs to
        T td = R_.cols[0].dot(moments[tgls[tid]].centroid);
        if ( td <= axdmp )  // add to 1st child
            tgls_[0].push_back(tgls[tid]);
        else
            tgls_[1].push_back(tgls[tid]);
    }

    const Vector3<T> cc = (minval + maxval)*0.5;
    c_ = R_ * cc;
    d_ = (maxval - minval) * 0.5;
}

/*
 * - eigen decomposition
 * - sort it based on eigenvalues
 */
template <typename T, class TTree>
void OBBTreeNode<T, TTree>::eigen_decomp(const Matrix3<T>& CM)
{
    T V[3][3], d[3];
    eigen_decomposition((const T(*)[3])&CM, V, d);

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
    R_.set(V[0][0], V[0][1], V[0][2],
            V[1][0], V[1][1], V[1][2],
            V[2][0], V[2][1], V[2][2]);
    assert((R_.cols[0].length_sqr()-1.) < 1E-9 &&
           (R_.cols[1].length_sqr()-1.) < 1E-9 &&
           (R_.cols[2].length_sqr()-1.) < 1E-9);
}

template <typename T, class TTree>
T OBBTreeNode<T, TTree>::dist_sqr_to_node(const Point3<T>& pt) const
{
    // coordinate in the local coordinate system of this node
    const Vector3<T> v = pt - c_;
    T ret = 0;

    Vector3<T> coord(R_.cols[0].dotProduct(v),
                     R_.cols[1].dotProduct(v),
                     R_.cols[2].dotProduct(v));
    for(int i = 0;i < 3;++ i)
    {
        if ( coord[i] < -d_[i] ) 
            ret += M_SQR(coord[i] + d_[i]);
        else if ( coord[i] >  d_[i] )
            ret += M_SQR(coord[i] - d_[i]);
    }
    return ret;
}

template <typename T, class TTree>
int OBBTreeNode<T, TTree>::is_disjoint(const OBBTreeNode<T, TTree>* nb) const
{
    const Vector3<T> ctrB2A = nb->c_ - c_;

    // check for a separation plane parallel to the face of A
    for(int ii = 0;ii < 3;++ ii)
    {
        // the projection direction is R_.cols[ii]
        const Vector3<T>& pd = R_.cols[ii];
        const T ctrProj = M_ABS(pd.dot(ctrB2A));
        const T dA = d_[ii];
        const T dB = M_ABS(pd.dot(nb->R_.cols[0]))*nb->d_[0] + 
                     M_ABS(pd.dot(nb->R_.cols[1]))*nb->d_[1] + 
                     M_ABS(pd.dot(nb->R_.cols[2]))*nb->d_[2];
        if ( ctrProj >= dA + dB ) return 1;
    }

    for(int ii = 0;ii < 3;++ ii)
    {
        // the projection direction is nb->R_.cols[ii]
        const Vector3<T>& pd = nb->R_.cols[ii];
        const T ctrProj = M_ABS(pd.dot(ctrB2A));
        const T dA = M_ABS(pd.dot(R_.cols[0]))*d_[0] +
                     M_ABS(pd.dot(R_.cols[1]))*d_[1] +
                     M_ABS(pd.dot(R_.cols[2]))*d_[2];
        const T dB = nb->d_[ii];
        if ( ctrProj >= dA + dB ) return 2;
    }

    for(int ii = 0;ii < 3;++ ii)
    for(int jj = 0;jj < 3;++ jj)
    {
        Vector3<T> pd = R_.cols[ii].cross(nb->R_.cols[jj]);
        T lpd = pd.normalize2();

        if ( lpd < PrecisionType<T>::EPS ) continue;

        const T ctrProj = M_ABS(pd.dot(ctrB2A));
        const T dA = M_ABS(pd.dot(R_.cols[0]))*d_[0] +
                     M_ABS(pd.dot(R_.cols[1]))*d_[1] +
                     M_ABS(pd.dot(R_.cols[2]))*d_[2];
        const T dB = M_ABS(pd.dot(nb->R_.cols[0]))*nb->d_[0] + 
                     M_ABS(pd.dot(nb->R_.cols[1]))*nb->d_[1] + 
                     M_ABS(pd.dot(nb->R_.cols[2]))*nb->d_[2];
        if ( ctrProj >= dA + dB ) return 3;
    }

    return 0;
}

#ifdef USE_NAMESPACE
}
#endif
#endif
