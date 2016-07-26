#ifndef OBB_TREE_HPP
#   define OBB_TREE_HPP

#include <vector>
#include "OBBTreeNode.hpp"
#include "Triangle.hpp"
#include "utils/math.hpp"
#include "generic/null_type.hpp"

#ifdef USE_NAMESPACE
namespace sploosh
{
#endif

/*!
 * Implement the oriented bounding box
 */
template <typename T, class _TMesh, class _TData=sploosh::NullType>
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
        {  return root_; }

        TNode* root()
        {  return root_; }

        const std::vector<moment>& triangle_moments() const
        {  return moments_; }

    private:
        void update_moments();

    private:
        TNode*                      root_;
        std::vector<moment>         moments_;
};

///////////////////////////////////////////////////////////////////////////////
template <typename T, class _TMesh, class _TData>
const int OBBTree<T, _TMesh, _TData>::LEAF_SIZE = 4;

template <typename T, class _TMesh, class _TData>
void OBBTree<T, _TMesh, _TData>::init(const TMesh* pmesh)
{
    if ( pmesh == NULL ) 
    {
        root_ = NULL;
        return;
    }

    const std::vector< Tuple3ui >&  indices = pmesh->surface_indices();
    const std::vector< Point3<T> >& vtx     = pmesh->vertices();
    std::vector<int> tgls(indices.size());
    moments_.resize(indices.size());

    const T TS = 1./3.;
    const T TS2 = 1./12.;

    //// compute moment for each triangle
    for(size_t tid = 0;tid < tgls.size();++ tid)
    {
        tgls[tid] = tid;

        const Point3<T>& v0 = vtx[indices[tid][0]]; 
        const Point3<T>& v1 = vtx[indices[tid][1]]; 
        const Point3<T>& v2 = vtx[indices[tid][2]];
        moments_[tid].area = Triangle<T>::area(v0, v1, v2);
        moments_[tid].centroid = (v0 + v1 + v2) * TS;

        if ( moments_[tid].area < 1E-18 )
        {
            for(int k = 0;k < 3;++ k)
            {
                moments_[tid].C.cols[k][k] = 
                    M_SQR(v0[k]) + M_SQR(v1[k]) + M_SQR(v2[k]);
                for(int j = k+1;j < 3;++ j)
                    moments_[tid].C.cols[k][j] = moments_[tid].C.cols[j][k] = 
                        v0[k]*v0[j] + v1[k]*v1[j] + v2[k]*v2[j];
            }
        }
        else
        {
            for(int k = 0;k < 3;++ k)
            {
                moments_[tid].C.cols[k][k] = moments_[tid].area * TS2 * 
                    (9.*M_SQR(moments_[tid].centroid[k]) + 
                     M_SQR(v0[k]) + M_SQR(v1[k]) + M_SQR(v2[k]));
                for(int j = k+1;j < 3;++ j)
                {
                    moments_[tid].C.cols[k][j] = moments_[tid].C.cols[j][k] = 
                        moments_[tid].area * TS2 * 
                        (9.*moments_[tid].centroid[k]*moments_[tid].centroid[j] +
                         v0[k]*v0[j] + v1[k]*v1[j] + v2[k]*v2[j]);
                }
            }
        }
    }  // end for 

    root_ = new OBBTreeNode<T, TSelf>(this, pmesh, tgls, 0);
}

template <typename T, class _TMesh, class _TData>
OBBTree<T, _TMesh, _TData>::~OBBTree()
{
    delete root_;
}

#ifdef USE_NAMESPACE
}
#endif
#endif
