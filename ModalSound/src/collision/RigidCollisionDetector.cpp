#include "RigidCollisionDetector.h"
#include "geometry/tritri.h"

using namespace std;

/*
 * Check the triangle pair intersection
 */
void RigidCollisionDetector::tri_tri_intersection_test(
        const std::pair<const TOBBNode*, const TOBBNode*>& pns,
        std::vector< std::pair<int, int> >& ret) const
{
    const vector< Point3d >& vtxA = meshA_->vertices();
    const vector< Point3d >& vtxB = meshB_->vertices();

    for(int ii = 0;ii < 2;++ ii)
    {
        const vector<int>& triA = pns.first->triangles(ii);
        if ( triA.empty() ) continue;

        for(int jj = 0;jj < 2;++ jj)
        {
            const vector<int>& triB = pns.second->triangles(jj);

            for(size_t iB = 0;iB < triB.size();++ iB)
            for(size_t iA = 0;iA < triA.size();++ iA)
            {
                const Tuple3ui& vIdA = meshA_->triangle_ids( triA[iA] );
                const Tuple3ui& vIdB = meshB_->triangle_ids( triB[iB] );

                if ( tri_tri_overlap_test_3d(
                        (const double*)(&(vtxA[vIdA.x])),
                        (const double*)(&(vtxA[vIdA.y])),
                        (const double*)(&(vtxA[vIdA.z])),
                        (const double*)(&(vtxB[vIdB.x])),
                        (const double*)(&(vtxB[vIdB.y])),
                        (const double*)(&(vtxB[vIdB.z]))) )
                    ret.push_back(make_pair(triA[iA], triB[iB]));
            }
        } // end for
    } // end for
}

int RigidCollisionDetector::detect(std::vector< std::pair<int, int> >& ret) 
{
    ret.clear();
    const TOBBNode* nnA = treeA_.root();
    const TOBBNode* nnB = treeB_.root();


    bvStack_.push(make_pair(nnA, nnB));

    while ( !bvStack_.empty() )
    {
        std::pair<const TOBBNode*, const TOBBNode*> pns = bvStack_.top();
        bvStack_.pop();

        // bounding box intersect each other
        if ( !pns.first->is_disjoint(pns.second) )
        {
            if ( pns.first->is_leaf() && pns.second->is_leaf() )
            {
                // triangle-triangle intersection test
                tri_tri_intersection_test(pns, ret);
            }
            else
            {
                if ( !pns.first->is_leaf() )
                {
                    if ( !pns.second->is_leaf() &&
                         pns.first->size() < pns.second->size() )
                    {
                        bvStack_.push(make_pair(pns.first, pns.second->left_child()));
                        bvStack_.push(make_pair(pns.first, pns.second->right_child()));
                    }
                    else
                    {
                        bvStack_.push(make_pair(pns.first->left_child(),  pns.second));
                        bvStack_.push(make_pair(pns.first->right_child(), pns.second));
                    }
                }
                else
                {
                    bvStack_.push(make_pair(pns.first, pns.second->left_child()));
                    bvStack_.push(make_pair(pns.first, pns.second->right_child()));
                }
            }
        } // end if
    }

    return (int)ret.size();
}
