#ifndef RIGID_COLLISION_DETECT_INC
#   define RIGID_COLLISION_DETECT_INC

#include <algorithm>
#include <vector>
#include <stack>
#include "geometry/TriangleMesh.hpp"
#include "geometry/OBBTree.hpp"

class RigidCollisionDetector
{
    public:
        typedef OBBTree< double, TriangleMesh<double> > TOBBTree;
        typedef TOBBTree::TNode                     TOBBNode;

        RigidCollisionDetector(const TriangleMesh<double>* mA, const TriangleMesh<double>* mB):
                meshA_(mA), meshB_(mB), treeA_(mA), treeB_(mB)
        { }

        int detect(std::vector< std::pair<int, int> >& ret);

    private:
        void tri_tri_intersection_test(
                const std::pair<const TOBBNode*, const TOBBNode*>& pns,
                std::vector< std::pair<int, int> >& ret) const;

    private:
        const TriangleMesh<double>*   meshA_;
        const TriangleMesh<double>*   meshB_;
        const TOBBTree                treeA_;
        const TOBBTree                treeB_;

        std::stack< std::pair<const TOBBNode*, const TOBBNode*> >   bvStack_;
};

#endif
