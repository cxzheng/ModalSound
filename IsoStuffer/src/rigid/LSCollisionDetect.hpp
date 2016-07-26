/******************************************************************************
 *  File: LSCollisionDetect.hpp
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
#ifndef RIGID_LEVEL_SET_COLLISION_PROCESSOR_H
#   define RIGID_LEVEL_SET_COLLISION_PROCESSOR_H

#include "config.h"
#include <set>
#include <assert.h>

#ifdef USE_HASH_MAP
#   include <ext/hash_map>
#else
#   include <tr1/unordered_map>
#endif
#include "CollisionRec.hpp"
#include "DistFunc.hpp"
#include "geometry/OBBTree.hpp"
#include "geometry/TriangleMesh.hpp"
#include "levelset/LevelSet3.hpp"
#include "utils/math.hpp"
#include "CollisionConstraint.hpp"

template <typename T>
struct TreeNodeData
{
    int         ts;
    Point3<T>   predc0;     // center of the tree node in predicted configuration
    Matrix3<T>  predR;      // three principle direction of this tree node in
                            // predicted configuration
                            
    TreeNodeData():ts(0) { }
};

class RigidBodySimulator;
class ContactGraph;

template <typename T, class TRigidBody>
class CollisionProcessor
{
    template <typename T1, class TCollProc> friend class DistFunc;
    friend class RigidBodySimulator;
    friend class ContactGraph;

    public:
        typedef TRigidBody                              TRBody;
        typedef CollisionProcessor<T, TRigidBody>       TSelf;
        typedef typename TRigidBody::TMesh              TMesh;
        typedef TriangleMesh<T>                         TSurf;
        typedef OBBTree< T, TSurf, TreeNodeData<T> >    TTree;
        typedef TreeNodeData<T>                         TTreeNodeData;
        typedef typename TTree::TNode                   TTreeNode;
        typedef LevelSet3<T, DistFunc<T, TSelf>, true>  TLevelSet;
        typedef DistFunc<T, TSelf>                      TDistFunc;

    public:
        CollisionProcessor(TRigidBody* rbody, TMesh* pmesh);

        ~CollisionProcessor();

        /*
         * Find bodyA's deepest vertex that is inside of bodyB.
         * NOTE: make sure m_collidedVtx has got updated before calling this method
         *
         * if no collision detected, return false, otherwise, outRec is filled
         * about that deepest vertex
         */
        bool deepest_penetrating_vertex(TRigidBody* bodyB, CollisionRec<T>& outRec);

        bool is_colliding(TRigidBody* bodyB);

        /*
         * find the current deepest collision with walls, ground etc
         */
        bool detect_collision_constraint(
                const std::vector< CollisionConstraint<T, TSelf>* >& constraints,
                CollisionRec<T>& outRec);

        /*
         * \param vid is the vertex id in surface mesh
         */
        inline void add_collision_vtx_candidate(int vid)
        {   m_collidedVtx.insert(vid); }

        void clean_collision_vtx_candidates()
        {   m_collidedVtx.clear();  }

        /*
         * retrieve the triangle info for collision. The triangle is specified by
         * the triangle ID
         */
        const TriangleForCollision<T>& triangle_info_for_collision(int tid) const
        {  
            assert(tid < m_tglCollInfo.size() && tid >= 0);
            return m_tglCollInfo[tid];
        }

        const Vector3<T>& edge_pseudo_normal(int a, int b)
        {
            assert(a != b);

            if ( a > b ) std::swap(a, b);
            return m_edgePseudoNml[a][b];
        }

        /*
         * update the predicted position of the given vertex with id=vid
         */
        const Point3<T>& predicted_surf_vtx_position(int vid)
        {
            if ( m_vtxTimestamp[vid] == m_predTimestamp )
                return m_predVtx[vid];

            m_vtxTimestamp[vid] = m_predTimestamp;
            m_predVtx[vid] = mp_rigidObj->predicted_position(m_surface.vertex(vid));
            return m_predVtx[vid];
        }

        void update_tree_node_state(TreeNodeData<T>* data, const TTreeNode* node)
        {
            //// node->c() and node->R() are all in its initial configuration
            data->ts = m_predTimestamp;
            data->predc0 = mp_rigidObj->predicted_position(node->c());
            const Matrix3<T>& R0 = node->R();
            data->predR.set(mp_rigidObj->predicted_normal(R0.cols[0]),
                            mp_rigidObj->predicted_normal(R0.cols[1]),
                            mp_rigidObj->predicted_normal(R0.cols[2]));
        }

        void inc_pred_timestamp() 
        {  ++ m_predTimestamp; }

        /* compute the lumped area for each vertex */
        inline void update_surface_vtx_areas()
        {   m_surface.update_vertex_areas(); }

        /* ============= retrieval methods ============ */
        const TTree* bounding_tree() const
        {   return mp_boundingTree; }
        const TRigidBody* rigid_body() const
        {   return mp_rigidObj; }
        const TriangleMesh<T>* surface_mesh() const
        {   return &m_surface; }
        TTreeNode* bounding_tree_root()
        {   return mp_boundingTree->root(); }
        const TTreeNode* bounding_tree_root() const
        {   return mp_boundingTree->root(); }
        int pred_timestamp() const
        {   return m_predTimestamp; }
        const std::vector< Vector3<T> >& vtx_pseudo_nml() const
        {   return m_vtxPseudoNml; }

    private:
        void compute_tgl_pseudo_normals();
        void compute_vtx_pseudo_normals();
        void compute_edge_pseudo_normals();
        void compute_tgl_collision_info();

    private:
        TRigidBody*                         mp_rigidObj;
        TMesh*                              mp_mesh;        // pointer to the tet mesh
        TriangleMesh<T>                     m_surface;      // surface triangle mesh

        /* 
         * OBB is created at initial configuration, which is in world coordinate
         * OBB tree is created based on the surface mesh
         */
        TTree*                              mp_boundingTree;
        TLevelSet*                          mp_levelset;

        /*
         * When we detect if this rigid obj is colliding with another one,
         * this set is to maintain a set of vertices of THIS obj, which might 
         * have penetrated into the other object. These vertices are just candidates
         * got from OBB-tree overlap test.
         *
         * The ID maintained in m_collidedVtx is the vertex ID for surface mesh
         */
        std::set<int>                       m_collidedVtx;

        /*
         * m_predVtx stores the predicted position of each vertex. Each predicted 
         * position is evaluated lazily. m_predTimestamp and m_vtxTimestamp are used
         * to indicate if the vertex is dirty. i.e. if m_predTimestamp > m_vtxTimestamp[i],
         * the predicted position of vertex i needs to be re-computed
         *
         * m_predTimestamp is to determine whether or not the m_predVtx and obbtreenode->data
         * need to be updated
         */
        int                                 m_predTimestamp;
        std::vector<int>                    m_vtxTimestamp;
        /*
         * The predicted positions of surface vertices
         */
        std::vector< Point3<T> >            m_predVtx;

        // the pre-computed pseudo normal information, and triangle transformation info for computing
        // the signed distance field. 
        // (Implement the paper: Generating Signed Distance Fields From Triangle Meshes)
        // NOTE that all the information are computed at the rigid obj's initial configuration.
        std::vector< Vector3<T> >               m_tglPseudoNml; // pseudo normal at each triangle
        std::vector< Vector3<T> >               m_vtxPseudoNml; // pseudo normal on each vertex
#ifdef USE_UNORDERED_MAP
        std::tr1::unordered_map< int, std::tr1::unordered_map< int, Vector3<T> > >  m_edgePseudoNml;
#else
        __gnu_cxx::hash_map< int, __gnu_cxx::hash_map< int, Vector3<T> > >          m_edgePseudoNml;
#endif
        std::vector< TriangleForCollision<T> >  m_tglCollInfo;  // collision information for each triangle
};

///////////////////////////////////////////////////////////////////////////////

template <typename T, class TRigidBody>
CollisionProcessor<T, TRigidBody>::CollisionProcessor(TRigidBody* rbody, TMesh* pmesh):
        mp_rigidObj(rbody), mp_mesh(pmesh), m_predTimestamp(1)
{
    //// extract triangle mesh for the surface
    pmesh->extract_surface(m_surface);
    mp_boundingTree = new TTree(&m_surface);

    // find the min/max point of mesh
    Point3<T> minPt(1E+100, 1E+100, 1E+100);
    Point3<T> maxPt(-1E+100, -1E+100, -1E+100);
    const std::vector< Point3<T> >& vtx = m_surface.vertices();

    for(size_t i = 0;i < vtx.size();++ i)
    {
        minPt.x = fmin(minPt.x, vtx[i].x);
        minPt.y = fmin(minPt.y, vtx[i].y);
        minPt.z = fmin(minPt.z, vtx[i].z);

        maxPt.x = fmax(maxPt.x, vtx[i].x);
        maxPt.y = fmax(maxPt.y, vtx[i].y);
        maxPt.z = fmax(maxPt.z, vtx[i].z);
    }
    mp_levelset = new TLevelSet(minPt, maxPt, DistFunc<T, TSelf>(this));

    //// allocate space
    m_vtxTimestamp.resize(m_surface.num_vertices());
    m_predVtx.resize(m_surface.num_vertices());

    //// pre-compute pseudo-normals
    compute_tgl_pseudo_normals();
    compute_edge_pseudo_normals();
    compute_vtx_pseudo_normals();
    compute_tgl_collision_info();
}

template <typename T, class TRigidBody>
void CollisionProcessor<T, TRigidBody>::compute_tgl_collision_info()
{
    const std::vector< Point3<T> >& vtx = m_surface.vertices(); 
    const std::vector< Tuple3ui >&  tgl = m_surface.surface_indices();

    m_tglCollInfo.resize(tgl.size());
    for(size_t i = 0;i < tgl.size();++ i)
    for(int j = 0;j < 3;++ j)
    {
        int vid0 = tgl[i][j], vid1 = tgl[i][(j+1)%3], vid2 = tgl[i][(j+2)%3];
        m_tglCollInfo[i].vtxId[j] = vid0;
        m_tglCollInfo[i].sidelen[j] = vtx[vid0].distance(vtx[vid1]);

        Vector3<T> vec01 = vtx[vid1] - vtx[vid0];
        m_tglCollInfo[i].QT[j].cols[0] = vec01 / m_tglCollInfo[i].sidelen[j];
        m_tglCollInfo[i].QT[j].cols[2] = vec01.crossProduct(vtx[vid2] - vtx[vid0]);
        m_tglCollInfo[i].QT[j].cols[2].normalize();
        m_tglCollInfo[i].QT[j].cols[1] = m_tglCollInfo[i].QT[j].cols[2].crossProduct(
                m_tglCollInfo[i].QT[j].cols[0]);

        m_tglCollInfo[i].Q[j] = m_tglCollInfo[i].QT[j].transpose();
        m_tglCollInfo[i].x0[j] = -(m_tglCollInfo[i].Q[j] * vtx[vid0]);
    }
}

template <typename T, class TRigidBody>
void CollisionProcessor<T, TRigidBody>::compute_tgl_pseudo_normals()
{
    const std::vector< Point3<T> >& vtx = m_surface.vertices(); 
    const std::vector< Tuple3ui >&  tgl = m_surface.surface_indices();

    m_tglPseudoNml.resize(tgl.size());
    for(size_t i = 0;i < tgl.size();++ i)
    {
        m_tglPseudoNml[i] = Triangle<T>::normal(
                vtx[tgl[i][0]], vtx[tgl[i][1]], vtx[tgl[i][2]]);
        if ( m_tglPseudoNml[i].lengthSqr() < 1E-24 )
        {
            fprintf(stderr, "ERROR: triangle has zero area: %.30lf\n",
                    m_tglPseudoNml[i].lengthSqr());
            exit(1);
        }
        m_tglPseudoNml[i].normalize();
    }
}

template <typename T, class TRigidBody>
void CollisionProcessor<T, TRigidBody>::compute_edge_pseudo_normals()
{
    m_edgePseudoNml.clear();

    const std::vector< Tuple3ui >&  tgl = m_surface.surface_indices();
    for(size_t i = 0;i < tgl.size();++ i)
    {
        int a, b;
        // edge 0
        a = std::min(tgl[i][0], tgl[i][1]);
        b = std::max(tgl[i][0], tgl[i][1]);
        m_edgePseudoNml[a][b] += m_tglPseudoNml[i];

        // edge 1
        a = std::min(tgl[i][1], tgl[i][2]);
        b = std::max(tgl[i][1], tgl[i][2]);
        m_edgePseudoNml[a][b] += m_tglPseudoNml[i];

        // edge 1
        a = std::min(tgl[i][2], tgl[i][0]);
        b = std::max(tgl[i][2], tgl[i][0]);
        m_edgePseudoNml[a][b] += m_tglPseudoNml[i];
    }
}

template <typename T, class TRigidBody>
void CollisionProcessor<T, TRigidBody>::compute_vtx_pseudo_normals()
{
    const std::vector< Point3<T> >& vtx = m_surface.vertices(); 
    const std::vector< Tuple3ui >&  tgl = m_surface.surface_indices();

    m_vtxPseudoNml.resize(vtx.size());
    memset(&m_vtxPseudoNml[0], 0, sizeof(Vector3<T>)*m_vtxPseudoNml.size());

    for(size_t i = 0;i < tgl.size();++ i)
    {
        const Vector3<T>& nml = m_tglPseudoNml[i];

        m_vtxPseudoNml[tgl[i][0]] += nml * Triangle<T>::angle(
                vtx[tgl[i][2]], vtx[tgl[i][0]], vtx[tgl[i][1]]);
        m_vtxPseudoNml[tgl[i][1]] += nml * Triangle<T>::angle(
                vtx[tgl[i][0]], vtx[tgl[i][1]], vtx[tgl[i][2]]);
        m_vtxPseudoNml[tgl[i][2]] += nml * Triangle<T>::angle(
                vtx[tgl[i][1]], vtx[tgl[i][2]], vtx[tgl[i][0]]);
    }
}

template <typename T, class TRigidBody>
CollisionProcessor<T, TRigidBody>::~CollisionProcessor()
{
    delete mp_boundingTree;
    delete mp_levelset;
}

/*
 * - go through each vertex recorded in m_collidedVtx
 * - ask bodyB->level_set() to compute the distance field value, and the normal
 *   direction.
 *   And find the deepest vertex that has a non-separating relative velocity
 *
 * NOTE: make sure m_collidedVtx has got updated before calling this method
 *
 * if no penetrating vertex is found, return false, otherwise, outRec is filled
 * about that deepest vertex
 */
template <typename T, class TRigidBody>
bool CollisionProcessor<T, TRigidBody>::deepest_penetrating_vertex(
        TRigidBody* bodyB, CollisionRec<T>& outRec)
{
    bool ret = false;
    Vector3<T> nml;
    outRec.depth = 1E+20;

    std::set<int>::iterator end = m_collidedVtx.end();
    for(std::set<int>::iterator it = m_collidedVtx.begin();
            it != end; ++ it)
    {
        // get the vertex position in predicted configuration, and
        // transform the vertex position from the other object's predicted
        // state into its initial state
        const Point3<T> predPt = predicted_surf_vtx_position(*it);
        const Point3<T> pt     = bodyB->initial_predicted_position(predPt);
        
        // check if the vertex has negative distance value
        const T isoval = bodyB->collision_processor()->mp_levelset->
                    negative_dist_with_normal(pt, nml);
        if ( isoval < 0. && isoval < outRec.depth )
        {
            // transform the nml into bodyB's current configuration
            const Vector3<T> predNml = bodyB->predicted_normal(nml);

            // check if the vertex and the object have nonseparating relative velocity
            // now the normal is point outward to the bodyB
            // (predPt - ptorigA) - (predPt - ptorigB) = ptorigB - ptorigA
            const Vector3<T> vab = mp_rigidObj->predicted_velocity(predPt) - 
                                    bodyB->predicted_velocity(predPt);
            const T vnrel = vab.dotProduct(predNml);

            if ( vnrel < 0. )   // nonseparating velocity
            {
                outRec.depth      = isoval;
                outRec.pt         = predPt;
                outRec.impulseDir = predNml;
                outRec.vnrel      = vnrel;
                outRec.vrel       = vab;
                outRec.eps        = fmin(bodyB->rest_coeff(),
                                         mp_rigidObj->rest_coeff());
                outRec.mu         = fmax(bodyB->friction_coeff(),
                                         mp_rigidObj->friction_coeff());
#ifdef USE_RECORDER
                outRec.vtxId      = *it;
#endif
                ret = true;
            } // end if
        }
    }

    return ret;
}

template <typename T, class TRigidBody>
bool CollisionProcessor<T, TRigidBody>::is_colliding(TRigidBody* bodyB)
{
    std::set<int>::iterator end = m_collidedVtx.end();
    for(std::set<int>::iterator it = m_collidedVtx.begin();
            it != end; ++ it)
    {
        // get the vertex position in predicted configuration, and
        // transform the vertex position from the other object's predicted
        // state into its initial state
        const Point3<T> predPt = predicted_surf_vtx_position(*it);
        const Point3<T> pt     = bodyB->initial_predicted_position(predPt);
        
        // check if the vertex has negative distance value
        if ( bodyB->collision_processor()->mp_levelset->distance(pt) < 0. ) 
            return true;
    }

    return false;
}

template <typename T, class TRigidBody>
bool CollisionProcessor<T, TRigidBody>::detect_collision_constraint(
        const std::vector< CollisionConstraint<T, TSelf>* >& constraints,
        CollisionRec<T>& outRec)
{

    bool ret = false;
    outRec.depth = 10.;
    for(size_t i = 0;i < constraints.size();++ i)
        ret |= constraints[i]->deepest_collision(this, 
                mp_boundingTree->root(), outRec);
    return ret;
}

///////////////////////////////////////////////////////////////////////////////////

/*
 * Do the overlap test as did in paper Gottschalk et.al. 1996
 * (OBBTree: A Hierarchical Structure for Rapid Interference Detection)
 */
template <typename T, class TRigidBody>
static bool is_disjoint(
        typename CollisionProcessor<T, TRigidBody>::TTreeNode* bdNodeA,
        CollisionProcessor<T, TRigidBody>* procA,
        typename CollisionProcessor<T, TRigidBody>::TTreeNode* bdNodeB,
        CollisionProcessor<T, TRigidBody>* procB)
{
    TreeNodeData<T>* dataA = bdNodeA->data();
    TreeNodeData<T>* dataB = bdNodeB->data();
    
    // make sure the predc0, predR are updated 
    if ( dataA->ts < procA->pred_timestamp() )
        procA->update_tree_node_state(dataA, bdNodeA);
    if ( dataB->ts < procB->pred_timestamp() )
        procB->update_tree_node_state(dataB, bdNodeB);

    const Vector3<T> dirab = dataA->predc0 - dataB->predc0;
    const Tuple3<T>& ra = bdNodeA->r();
    const Tuple3<T>& rb = bdNodeB->r();

    // the three principle dir of A as the projection dir
    for(int i = 0;i < 3;++ i)
    {
        T tl = fabs(dirab.dotProduct(dataA->predR.cols[i]));
        T tr = fabs(dataB->predR.cols[0].dotProduct(dataA->predR.cols[i])*rb[0]) + 
               fabs(dataB->predR.cols[1].dotProduct(dataA->predR.cols[i])*rb[1]) +
               fabs(dataB->predR.cols[2].dotProduct(dataA->predR.cols[i])*rb[2]) +
               ra[i];
        if ( tl > tr ) return true;
    }

    // the three principle dir of B as the projection dir
    for(int i = 0;i < 3;++ i)
    {
        T tl = fabs(dirab.dotProduct(dataB->predR.cols[i]));
        T tr = fabs(dataA->predR.cols[0].dotProduct(dataB->predR.cols[i])*ra[0]) +
               fabs(dataA->predR.cols[1].dotProduct(dataB->predR.cols[i])*ra[1]) +
               fabs(dataA->predR.cols[2].dotProduct(dataB->predR.cols[i])*ra[2]) +
               rb[i];
        if ( tl > tr ) return true;
    }

    double R[3][3];
    for(int i = 0;i < 3;++ i)
    for(int j = 0;j < 3;++ j)
        R[i][j] = fabs(dataA->predR.cols[i].dotProduct(dataB->predR.cols[j]));

    // cross product of the principle dir A and B as the projection dir
    for(int i = 0;i < 3;++ i)
    for(int j = 0;j < 3;++ j)
    {
        // pd = A_i x B_j
        Vector3<T> pd = dataA->predR.cols[i].crossProduct(dataB->predR.cols[j]);
        if ( pd.lengthSqr() < 1E-18 ) continue;

        T tl = fabs(dirab.dotProduct(pd));
        T tr = 0;
        for(int k = 0;k < 3;++ k)
        {
            if ( k != i ) tr += ra[k]*R[3-i-k][j];
            if ( k != j ) tr += rb[k]*R[i][3-j-k];
        }
        if ( tl > tr ) return true;
    }
    return false;
}

/*
 * Add all the vertices in nodeA as candidates colliding with bodyB
 * 
 * It also estimates a distance from each vertex in nodeA to bodyB by
 * check the minimum value
 *
 * both bdNodeA and bdNodeB are leaf node
 */
template <typename T, class TRigidBody>
static void collision_vtx_candiates(
        typename CollisionProcessor<T, TRigidBody>::TTreeNode* bdNodeA,
        CollisionProcessor<T, TRigidBody>* procA)
{
    /* vertex id for each triangle */
    const std::vector<Tuple3ui>&  tvIdsA = procA->surface_mesh()->surface_indices();

    for(int cid = 0;cid < 2;++ cid) // iterate on two sets of triangles
    {
        const std::vector<int>& ts = bdNodeA->triangles(cid);
        for(size_t i = 0;i < ts.size();++ i)    // for each triangle
        for(int vid = 0;vid < 3;++ vid)         // for each vertex of a triangle
            procA->add_collision_vtx_candidate(tvIdsA[ts[i]][vid]);
    }
}

/*
 * detect collisions of rigid objects
 *
 * collision detection is processed by
 * - using the precomputed OBB tree to find the leaf nodes where overlapping occurs
 * - go through all the vertices contained in each pair of OBB leaf nodes check if
 *   a vertex is penetrating into the other obj. This step is processed by looking
 *   up a discrete distance field table that is built lazily.
 */
template <typename T, class TRigidBody>
static void detect_tree_node_collisions(
        typename CollisionProcessor<T, TRigidBody>::TTreeNode* bdNodeA,
        typename CollisionProcessor<T, TRigidBody>::TTreeNode* bdNodeB,
        CollisionProcessor<T, TRigidBody>* procA,
        CollisionProcessor<T, TRigidBody>* procB)
{
    //// check if the bounding box in prediction state is overlap with each other
    if ( is_disjoint(bdNodeA, procA, bdNodeB, procB) ) return;
    
    //// if both bounding tree nodes are leaf nodes,
    //   do the level-set based collision detection
    if ( bdNodeA->is_leaf() && bdNodeB->is_leaf() )
    {
        collision_vtx_candiates(bdNodeA, procA);
        collision_vtx_candiates(bdNodeB, procB);
        return;
    }

    if ( !bdNodeA->is_leaf() )
    {
        // both A and B are not leaf
        //// divide the larger body
        if ( !bdNodeB->is_leaf() && bdNodeA->size() < bdNodeB->size() )
        {
            detect_tree_node_collisions(bdNodeA, 
                    (typename CollisionProcessor<T, TRigidBody>::TTreeNode *)bdNodeB->left_child(), 
                    procA, procB);
            detect_tree_node_collisions(bdNodeA, 
                    (typename CollisionProcessor<T, TRigidBody>::TTreeNode *)bdNodeB->right_child(), 
                    procA, procB);
        }
        else
        {
            detect_tree_node_collisions(
                    (typename CollisionProcessor<T, TRigidBody>::TTreeNode *)bdNodeA->left_child(), 
                    bdNodeB, procA, procB);
            detect_tree_node_collisions(
                    (typename CollisionProcessor<T, TRigidBody>::TTreeNode *)bdNodeA->right_child(),
                    bdNodeB, procA, procB);
        }
    }
    else    // node_a is leaf, node_b should not be the leaf
    {
        detect_tree_node_collisions(bdNodeA, 
                (typename CollisionProcessor<T, TRigidBody>::TTreeNode *)bdNodeB->left_child(), 
                procA, procB);
        detect_tree_node_collisions(bdNodeA, 
                (typename CollisionProcessor<T, TRigidBody>::TTreeNode *)bdNodeB->right_child(), 
                procA, procB);
    }
}

#endif

