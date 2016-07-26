#ifndef GEOMETRY_TETMESH_HPP
#   define GEOMETRY_TETMESH_HPP

#include "config.h"
#include <vector>
#include <valarray>
#include <set>
#include <boost/unordered_map.hpp>
#include <string.h>
#include "Tet.hpp"
#include "TriangleMesh.hpp"
#include "utils/tuple.hpp"
#include "sc/Tuple4.hpp"

#ifdef USE_NAMESPACE
namespace sploosh
{
#endif

/*
 * the velocity and acceleration for each individual vertex is enabled when the
 * macro WITH_DEFORMABLE_TET is defined
 */
template <typename T>
class TetMesh
{
    public:
        typedef Tuple4ui    TetIdx;  // vertices' index for each tet vertex

        // neighbor record for each tet
        struct NeighborRec
        {
            int num;                // the # of neighbor tets
            int id[4];              // the id of neighbor tets
        };

        /* copy constructor */
        TetMesh(const TetMesh<T>& rhs);

#if defined(DEBUG) | defined(_DEBUG)
        TetMesh():canAddVtx_(true) { }
#else
        TetMesh() { }
#endif

        /*
         * clone the given mesh without initialization
         */
        void uninitialized_clone(const TetMesh<T>* rhs);

        void clear();
        /*! add a tet with the indices of its 4 vertices */
        void add_tet(int v0, int v1, int v2, int v3);

        void update_surface();

        void extract_surface(TriangleMesh<T>* smesh) const; 
        void surface_id_map(std::map<int, int>& idmap) const;

        void update_tet_neighbors();
        void update_tet_neighbors(std::vector<NeighborRec>&) const;
        void update_tet_neighbors(
                const std::vector<bool>& tetOK,
                std::vector<NeighborRec>& neighbors) const;
        ///*
        // * get the tet list associated with each vertex
        // */
        //void get_vtx_tets_table(std::vector< std::vector<int> >& tbl) const;

        void get_vtx_neighbor_table(std::vector< std::set<int> >& tbl) const;

        void init();

        void reset_vertices()
        {
            memcpy(&vertices_[0], &restPos_[0], sizeof(Point3<T>)*restPos_.size());
        }

        template <typename FromT>
        void add_vertex(const Point3<FromT>& v)
        {
#if defined(DEBUG) | defined(_DEBUG)
            assert(canAddVtx_);
#endif

            vertices_.push_back(v);
            restPos_.push_back(v);
        }

        /*
         * add a vertex with returning the id of the newly added vertex
         */
        template <typename FromT>
        int add_vertex_rt(const Point3<FromT>& v)
        {
#if defined(DEBUG) | defined(_DEBUG)
            assert(canAddVtx_);
#endif

            vertices_.push_back(v);
            restPos_.push_back(v);

            return vertices_.size() - 1;
        }

        /*
         * ### NOTE: ###
         * Be aware of what you are doing when calling this method. It may 
         * screw up the pointers in tets_
         *
         * return the ID of the added vertex
         */
        template <typename FromT>
        int add_vertex_unsafe(const Point3<FromT>& now, const Point3<FromT>& rest)
        {
            restPos_.push_back(rest);
            vertices_.push_back(now);
            return vertices_.size() - 1;
        }

        void update_deformation_gradient();

        /*
         * If any new vertex is added after some tet has been added. The pointer in
         * vertices_ may point to wrong address. The method should be called to make
         * sure the pointer in Tet<T> maintains the correct address.
         *
         * void conform_vtx_pointers();
         */

        /* ============ retrival methods ============ */

        /*! return number of tetrahedrons */
        size_t num_tets() const { return tetIdx_.size(); }
        size_t num_vertices() const { return vertices_.size(); }
        size_t num_surface_tgls() const { return surfIdx_.size(); }

        const std::vector<TetIdx>& tet_indices() const { return tetIdx_; }

        const std::vector<Tuple3ui>& surface_indices() const
        { return surfIdx_; }

        /*! return masses associated with each vertices */
        const std::vector<T>& masses() const 
        { return masses_; }

        const std::vector< Tet<T> >& tets() const 
        { return tets_; }

        std::vector< Tet<T> >& tets() 
        { return tets_; }

        std::vector< Point3<T> >& vertices()
        { return vertices_; }

        const Point3<T>& vertex(size_t vid) const
        { 
            assert(vid < vertices_.size());
            return vertices_[vid]; 
        }

        const Point3<T>& rest_vertex(size_t vid) const
        {
            assert(vid < restPos_.size());
            return restPos_[vid]; 
        }

        const std::vector< Vector3<T> >& normals() const
        { return normals_; }

        const std::vector< Point3<T> >& vertices() const
        { return vertices_; }

        const std::vector< Point3<T> >& rest_positions() const
        { return restPos_; }

        std::vector< Point3<T> >& rest_positions()
        { return restPos_; }

#ifdef WITH_DEFORMABLE_TET
        std::vector< Vector3<T> >& velocities()
        { return vel_; }

        std::vector< Vector3<T> >& accelerations()
        { return acc_; }
#endif

        /*
         * surf2tet_ is computed when update_surface is called
         */
        const std::vector<int>& surf_to_tet() const
        { return surf2tet_; }

        /* 
         * check if the mesh is write referenced
         * only for debug purpose
         */
        bool check_mesh() const
        {
            for(size_t i = 0;i < tetIdx_.size();++ i)
            {
                for(int j = 0;j < 4;++ j)
                    if ( &(vertices_[tetIdx_[i][j]]) != tets_[i].vtx(j) ) 
                    {
                        fprintf(stderr, "WARNING: CHECK Mesh Failed at tet [%d] %d\n", 
                                (int)i, (int)tetIdx_[i][j]);
                        return false;
                    }
            }
            return true;
        }

    private:
        void add_normal(unsigned int a, unsigned int b, unsigned int c,
                std::vector<int>& cnts);

    protected:
        std::vector< Point3<T> >    vertices_;     // current positions of vertices
        std::vector< Point3<T> >    restPos_;      // position at rest pose

    private:
#if defined(DEBUG) | defined(_DEBUG)
        bool        canAddVtx_;
#endif

#ifdef WITH_DEFORMABLE_TET
        std::vector< Vector3<T> >   vel_;          // velocity of each vertex
        std::vector< Vector3<T> >   acc_;          // acceleration of each vertex
#endif
        std::vector< Tet<T> >       tets_;
        std::vector< TetIdx >       tetIdx_;       // vertex index for each vertex of the tet
        std::vector< Tuple3ui >     surfIdx_;      // vertices' index for surface triangles
        /*!
         * the vector<int> which maps each triangle in surfIdx_ to
         * tet that has the triangle
         */
        std::vector<int>            surf2tet_;

        std::vector<T>              masses_;       // masses_[i] is the mass of vertices_[i]
        std::vector< Vector3<T> >   normals_;
};

///////////////////////////////////////////////////////////////////////////////

template <typename T>
TetMesh<T>::TetMesh(const TetMesh<T>& rhs)
{
    clear();
    for(size_t i = 0;i < rhs.vertices_.size();++ i)
        add_vertex_unsafe(rhs.vertices_[i], rhs.restPos_[i]);

    for(size_t i = 0;i < rhs.tets_.size();++ i)
        add_tet(rhs.tetIdx_[i][0], rhs.tetIdx_[i][1], 
                rhs.tetIdx_[i][2], rhs.tetIdx_[i][3]);

#ifdef WITH_DEFORMABLE_TET
    //// initialize it
    vel_.resize(vertices_.size());
    acc_.resize(vertices_.size());
    memset(&vel_[0], 0, sizeof(Vector3<T>)*vel_.size());
    memset(&acc_[0], 0, sizeof(Vector3<T>)*acc_.size());
#endif

    masses_.resize(vertices_.size());
    tets_.resize(tetIdx_.size());
    for(size_t i = 0;i < tetIdx_.size();++ i)
    {
        tets_[i].init(&restPos_[tetIdx_[i][0]], &restPos_[tetIdx_[i][1]],
                       &restPos_[tetIdx_[i][2]], &restPos_[tetIdx_[i][3]]);
        T v = tets_[i].volume() * 0.25;

        tets_[i].conform_vtx_pointers(
                &vertices_[tetIdx_[i][0]],
                &vertices_[tetIdx_[i][1]],
                &vertices_[tetIdx_[i][2]],
                &vertices_[tetIdx_[i][3]]);

        masses_[tetIdx_[i][0]] += v;
        masses_[tetIdx_[i][1]] += v;
        masses_[tetIdx_[i][2]] += v;
        masses_[tetIdx_[i][3]] += v;
    }
#if defined(DEBUG) | defined(_DEBUG)
    canAddVtx_ = false;
#endif
}

template <typename T>
void TetMesh<T>::uninitialized_clone(const TetMesh<T>* rhs)
{
    clear();

    vertices_.resize(rhs->vertices_.size());
    restPos_.resize(rhs->restPos_.size());
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(rhs)
#endif
    for(int i = 0;i < (int)rhs->vertices_.size();++ i)
    {
        vertices_[i] = rhs->vertices_[i];
        restPos_[i]  = rhs->restPos_[i];
    }

    tetIdx_.resize(rhs->tetIdx_.size());
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(rhs)
#endif
    for(int i = 0;i < (int)rhs->tetIdx_.size();++ i)
        tetIdx_[i] = rhs->tetIdx_[i];

#ifdef WITH_DEFORMABLE_TET
    //// initialize it
    vel_.resize(vertices_.size());
    acc_.resize(vertices_.size());
    memset(&vel_[0], 0, sizeof(Vector3<T>)*vel_.size());
    memset(&acc_[0], 0, sizeof(Vector3<T>)*acc_.size());
#endif

    surfIdx_.resize(rhs->surfIdx_.size());
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(rhs)
#endif
    for(int i = 0;i < (int)rhs->surfIdx_.size();++ i)
        surfIdx_[i] = rhs->surfIdx_[i];

    if ( !rhs->surf2tet_.empty() )
    {
        surf2tet_.resize(rhs->surf2tet_.size());
#ifdef USE_OPENMP
        #pragma omp parallel for default(none) shared(rhs)
#endif
        for(int i = 0;i < (int)rhs->surf2tet_.size();++ i)
            surf2tet_[i] = rhs->surf2tet_[i];
    }

    if ( !rhs->normals_.empty() )
    {
        normals_.resize(rhs->normals_.size());
#ifdef USE_OPENMP
        #pragma omp parallel for default(none) shared(rhs)
#endif
        for(int i = 0;i < (int)rhs->normals_.size();++ i)
            normals_[i] = rhs->normals_[i];
    }
}

/*!
 * - initialize the velocity and acceleration vector
 * - compute the mass and inverse mass vector
 *
 * NOTE: we create the tets_ array here, assuming that no add_vertex is gonna call
 *       after this. Otherwise, a growing m_vertex array may change its address, and 
 *       hence screw up the vertex pointers in m_tet
 */
template <typename T>
void TetMesh<T>::init()
{
#ifdef WITH_DEFORMABLE_TET
    vel_.resize(vertices_.size());
    acc_.resize(vertices_.size());
    memset(&vel_[0], 0, sizeof(Vector3<T>)*vel_.size());
    memset(&acc_[0], 0, sizeof(Vector3<T>)*acc_.size());
#endif

    masses_.resize(vertices_.size());
    tets_.resize(tetIdx_.size());
    for(size_t i = 0;i < tetIdx_.size();++ i)
    {
        tets_[i].init(&vertices_[tetIdx_[i][0]],
                &vertices_[tetIdx_[i][1]],
                &vertices_[tetIdx_[i][2]],
                &vertices_[tetIdx_[i][3]]);
        T v = tets_[i].volume() * 0.25;

        masses_[tetIdx_[i][0]] += v;
        masses_[tetIdx_[i][1]] += v;
        masses_[tetIdx_[i][2]] += v;
        masses_[tetIdx_[i][3]] += v;
    }
#if defined(DEBUG) | defined(_DEBUG)
    canAddVtx_ = false;
#endif
}

/*
 * NOTE: don't call 
 *  tets_.push_back(Tet<T>(&vertices_[v0], &vertices_[v1],
 *              &vertices_[v2], &vertices_[v3]));
 * here. This is because vertices_ is still growing. System might 
 * change the address of vertices_
 */
template <typename T>
void TetMesh<T>::add_tet(int v0, int v1, int v2, int v3)
{
    assert(v0 < (int)vertices_.size() && v1 < (int)vertices_.size() &&
           v2 < (int)vertices_.size() && v3 < (int)vertices_.size());

    tetIdx_.push_back(TetIdx(v0, v1, v2, v3));
}

template <typename T>
void TetMesh<T>::clear()
{
    vertices_.clear();
    restPos_.clear();
#ifdef WITH_DEFORMABLE_TET    
    vel_.clear();
    acc_.clear();
#endif
    tets_.clear();
    tetIdx_.clear();
    surfIdx_.clear();
    masses_.clear();
    normals_.clear();
}

/*!
 * NOTE: if size of tetRemoved is less than the num of tets, only the
 *       first "tetRemoved.size()" tets in this tet mesh is labeled by the
 *       tetRemoved array, all the later tet are assumed NOT to be removed by
 *       default.
 */
template <typename T>
void TetMesh<T>::update_tet_neighbors(
        const std::vector<bool>& tetRemoved,
        std::vector<NeighborRec>& neighbors) const
{
    using namespace boost;

    int v0, v1, v2;
    unordered_map< int, unordered_map<int, unordered_map<int, int> > > hash;

    neighbors.resize(tetIdx_.size());
    for(size_t i = 0;i < tetIdx_.size(); ++ i)
        neighbors[i].num = 0;

    for(size_t i = 0;i < tetIdx_.size();++ i)
    {
        if ( i < tetRemoved.size() && tetRemoved[i] ) continue;

        // face 0
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);

        if ( hash.count(v0) && hash[v0].count(v1) && hash[v0][v1].count(v2) )
        {
            int nid = hash[v0][v1][v2];
            neighbors[i].id[neighbors[i].num ++] = nid;
            neighbors[nid].id[neighbors[nid].num ++] = i;
        }
        else
            hash[v0][v1][v2] = i;

        // face 1
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][2];
        sort_triple(v0, v1, v2);

        if ( hash.count(v0) && hash[v0].count(v1) && hash[v0][v1].count(v2) )
        {
            int nid = hash[v0][v1][v2];
            neighbors[i].id[neighbors[i].num ++] = nid;
            neighbors[nid].id[neighbors[nid].num ++] = i;
        }
        else
            hash[v0][v1][v2] = i;

        // face 2
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);

        if ( hash.count(v0) && hash[v0].count(v1) && hash[v0][v1].count(v2) )
        {
            int nid = hash[v0][v1][v2];
            neighbors[i].id[neighbors[i].num ++] = nid;
            neighbors[nid].id[neighbors[nid].num ++] = i;
        }
        else
            hash[v0][v1][v2] = i;

        // face 3
        v0 = tetIdx_[i][1];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);

        if ( hash.count(v0) && hash[v0].count(v1) && hash[v0][v1].count(v2) )
        {
            int nid = hash[v0][v1][v2];
            neighbors[i].id[neighbors[i].num ++] = nid;
            neighbors[nid].id[neighbors[nid].num ++] = i;
        }
        else
            hash[v0][v1][v2] = i;
    }

}

template <typename T>
void TetMesh<T>::update_tet_neighbors(std::vector<NeighborRec>& neighbors) const
{
    using namespace boost;

    int v0, v1, v2;
    unordered_map< int, unordered_map<int, unordered_map<int, int> > > hash;

    neighbors.resize(tetIdx_.size());
    for(size_t i = 0;i < tetIdx_.size(); ++ i)
        neighbors[i].num = 0;

    for(size_t i = 0;i < tetIdx_.size();++ i)
    {
        // face 0
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash.count(v0) && hash[v0].count(v1) && hash[v0][v1].count(v2) )
        {
            int nid = hash[v0][v1][v2];
            neighbors[i].id[neighbors[i].num ++] = nid;
            neighbors[nid].id[neighbors[nid].num ++] = i;
        }
        else
            hash[v0][v1][v2] = i;

        // face 1
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][2];
        sort_triple(v0, v1, v2);
        if ( hash.count(v0) && hash[v0].count(v1) && hash[v0][v1].count(v2) )
        {
            int nid = hash[v0][v1][v2];
            neighbors[i].id[neighbors[i].num ++] = nid;
            neighbors[nid].id[neighbors[nid].num ++] = i;
        }
        else
            hash[v0][v1][v2] = i;

        // face 2
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash.count(v0) && hash[v0].count(v1) && hash[v0][v1].count(v2) )
        {
            int nid = hash[v0][v1][v2];
            neighbors[i].id[neighbors[i].num ++] = nid;
            neighbors[nid].id[neighbors[nid].num ++] = i;
        }
        else
            hash[v0][v1][v2] = i;

        // face 3
        v0 = tetIdx_[i][1];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash.count(v0) && hash[v0].count(v1) && hash[v0][v1].count(v2) )
        {
            int nid = hash[v0][v1][v2];
            neighbors[i].id[neighbors[i].num ++] = nid;
            neighbors[nid].id[neighbors[nid].num ++] = i;
        }
        else
            hash[v0][v1][v2] = i;
    }
}

/*!
 * Update the surface of the tet mesh. O(n.log(n)) complexity
 * and compute the normal for each surface vertex
 */
template <typename T>
void TetMesh<T>::update_surface()
{
    using namespace boost;

    int v0, v1, v2;
    unordered_map< int, unordered_map<int, unordered_map<int, int> > > hash;
    unordered_map< int, unordered_map<int, unordered_map<int, int> > > hash2;

    if ( normals_.size() != vertices_.size() )
        normals_.resize(vertices_.size());

#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 2000)
    for(size_t i = 0;i < vertices_.size();++ i)
#else
    for(int i = 0;i < (int)vertices_.size();++ i)
#endif
        normals_[i].zero();

    for(size_t i = 0;i < tetIdx_.size();++ i)
    {
        // face 0
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
        hash2[v0][v1][v2] = i;

        // face 1
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][2];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
        hash2[v0][v1][v2] = i;

        // face 2
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
        hash2[v0][v1][v2] = i;

        // face 3
        v0 = tetIdx_[i][1];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
        hash2[v0][v1][v2] = i;
    }

    surfIdx_.clear();
    surf2tet_.clear();

    std::vector<int> cnts(vertices_.size());
    memset(&cnts[0], 0, sizeof(int)*cnts.size());

    for(size_t i = 0;i < tetIdx_.size();++ i)
    {
        // face 0
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            surfIdx_.push_back(Tuple3i(
                        tetIdx_[i][0], 
                        tetIdx_[i][1],
                        tetIdx_[i][3]));
            surf2tet_.push_back(hash2[v0][v1][v2]);
            add_normal(tetIdx_[i][0], tetIdx_[i][1], tetIdx_[i][3], cnts);
        }

        // face 1
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][2];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            surfIdx_.push_back(Tuple3i(
                        tetIdx_[i][0], 
                        tetIdx_[i][2],
                        tetIdx_[i][1]));
            surf2tet_.push_back(hash2[v0][v1][v2]);
            add_normal(tetIdx_[i][0], tetIdx_[i][2], tetIdx_[i][1], cnts);
        }

        // face 2
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            surfIdx_.push_back(Tuple3i(
                        tetIdx_[i][3], 
                        tetIdx_[i][2],
                        tetIdx_[i][0]));
            surf2tet_.push_back(hash2[v0][v1][v2]);
            add_normal(tetIdx_[i][3], tetIdx_[i][2], tetIdx_[i][0], cnts);
        }

        // face 3
        v0 = tetIdx_[i][1];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            surfIdx_.push_back(Tuple3i(
                        tetIdx_[i][1], 
                        tetIdx_[i][2],
                        tetIdx_[i][3]));
            surf2tet_.push_back(hash2[v0][v1][v2]);
            add_normal(tetIdx_[i][1], tetIdx_[i][2], tetIdx_[i][3], cnts);
        }
    }

#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 2000) shared(cnts)
    for(size_t i = 0;i < vertices_.size();++ i)
#else
    for(int i = 0;i < (int)vertices_.size();++ i)
#endif
    {
        normals_[i] /= (T)cnts[i];
        normals_[i].normalize();
    }
    printf("Update surface: %d triangles\n", (int)surfIdx_.size());
}

template <typename T>
void TetMesh<T>::update_deformation_gradient()
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 2000)
    for(int i = 0;i < (int)tets_.size();++ i)
#else
    for(size_t i = 0;i < tets_.size();++ i)
#endif
        tets_[i].deformation_gradient();
}

template <typename T>
void TetMesh<T>::add_normal(unsigned int a, unsigned int b, 
        unsigned int c, std::vector<int>& cnts)
{
    Vector3<T> n = (restPos_[b] - restPos_[a]).cross(
                restPos_[c] - restPos_[a]);

    ++ cnts[a];
    ++ cnts[b];
    ++ cnts[c];
    normals_[a] += n;
    normals_[b] += n;
    normals_[c] += n;
}

template <typename T>
void TetMesh<T>::extract_surface(TriangleMesh<T>* smesh) const
{
    using namespace boost;

    unordered_map< int, unordered_map<int, unordered_map<int, int> > > hash;
    unordered_map< int, int > idmap;

    smesh->clear();

    int v0, v1, v2;
    for(size_t i = 0;i < tetIdx_.size();++ i)
    {
        // face 0
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 1
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][2];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 2
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 3
        v0 = tetIdx_[i][1];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
    }

    for(size_t i = 0;i < tetIdx_.size();++ i)
    {
        // face 0
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(tetIdx_[i][0]) )
                idmap[tetIdx_[i][0]] = smesh->add_vertex(restPos_[tetIdx_[i][0]]);
            if ( !idmap.count(tetIdx_[i][1]) )
                idmap[tetIdx_[i][1]] = smesh->add_vertex(restPos_[tetIdx_[i][1]]);
            if ( !idmap.count(tetIdx_[i][3]) )
                idmap[tetIdx_[i][3]] = smesh->add_vertex(restPos_[tetIdx_[i][3]]);
            smesh->add_triangle(idmap[tetIdx_[i][0]], idmap[tetIdx_[i][1]], idmap[tetIdx_[i][3]]);
        }

        // face 1
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][2];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(tetIdx_[i][0]) )
                idmap[tetIdx_[i][0]] = smesh->add_vertex(restPos_[tetIdx_[i][0]]);
            if ( !idmap.count(tetIdx_[i][2]) )
                idmap[tetIdx_[i][2]] = smesh->add_vertex(restPos_[tetIdx_[i][2]]);
            if ( !idmap.count(tetIdx_[i][1]) )
                idmap[tetIdx_[i][1]] = smesh->add_vertex(restPos_[tetIdx_[i][1]]);
            smesh->add_triangle(idmap[tetIdx_[i][0]], idmap[tetIdx_[i][2]], idmap[tetIdx_[i][1]]);
        }

        // face 2
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(tetIdx_[i][0]) )
                idmap[tetIdx_[i][0]] = smesh->add_vertex(restPos_[tetIdx_[i][0]]);
            if ( !idmap.count(tetIdx_[i][3]) )
                idmap[tetIdx_[i][3]] = smesh->add_vertex(restPos_[tetIdx_[i][3]]);
            if ( !idmap.count(tetIdx_[i][2]) )
                idmap[tetIdx_[i][2]] = smesh->add_vertex(restPos_[tetIdx_[i][2]]);
            smesh->add_triangle(idmap[tetIdx_[i][0]], idmap[tetIdx_[i][3]], idmap[tetIdx_[i][2]]);
        }

        // face 3
        v0 = tetIdx_[i][1];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(tetIdx_[i][1]) )
                idmap[tetIdx_[i][1]] = smesh->add_vertex(restPos_[tetIdx_[i][1]]);
            if ( !idmap.count(tetIdx_[i][2]) )
                idmap[tetIdx_[i][2]] = smesh->add_vertex(restPos_[tetIdx_[i][2]]);
            if ( !idmap.count(tetIdx_[i][3]) )
                idmap[tetIdx_[i][3]] = smesh->add_vertex(restPos_[tetIdx_[i][3]]);
            smesh->add_triangle(idmap[tetIdx_[i][1]], idmap[tetIdx_[i][2]], idmap[tetIdx_[i][3]]);
        }
    }
}

/* 
 * Extract surface and get the idmap which maps the vertex id from tet mesh 
 * to triangle mesh
 */
template <typename T>
void TetMesh<T>::surface_id_map(std::map<int, int>& idmap) const
{
    using namespace boost;

    int v0, v1, v2;
    unordered_map< int, unordered_map<int, unordered_map<int, int> > > hash;

    idmap.clear();
    for(size_t i = 0;i < tetIdx_.size();++ i)
    {
        // face 0
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 1
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][2];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 2
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 3
        v0 = tetIdx_[i][1];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
    }

    int vtxCnt = 0;
    for(size_t i = 0;i < tetIdx_.size();++ i)
    {
        // face 0
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(tetIdx_[i][0]) )
                idmap[tetIdx_[i][0]] = vtxCnt ++;
            if ( !idmap.count(tetIdx_[i][1]) )
                idmap[tetIdx_[i][1]] = vtxCnt ++;
            if ( !idmap.count(tetIdx_[i][3]) )
                idmap[tetIdx_[i][3]] = vtxCnt ++;
        }

        // face 1
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][1];
        v2 = tetIdx_[i][2];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(tetIdx_[i][0]) )
                idmap[tetIdx_[i][0]] = vtxCnt ++;
            if ( !idmap.count(tetIdx_[i][2]) )
                idmap[tetIdx_[i][2]] = vtxCnt ++;
            if ( !idmap.count(tetIdx_[i][1]) )
                idmap[tetIdx_[i][1]] = vtxCnt ++;
        }

        // face 2
        v0 = tetIdx_[i][0];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(tetIdx_[i][0]) )
                idmap[tetIdx_[i][0]] = vtxCnt ++;
            if ( !idmap.count(tetIdx_[i][3]) )
                idmap[tetIdx_[i][3]] = vtxCnt ++;
            if ( !idmap.count(tetIdx_[i][2]) )
                idmap[tetIdx_[i][2]] = vtxCnt ++;
        }

        // face 3
        v0 = tetIdx_[i][1];
        v1 = tetIdx_[i][2];
        v2 = tetIdx_[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(tetIdx_[i][1]) )
                idmap[tetIdx_[i][1]] = vtxCnt ++;
            if ( !idmap.count(tetIdx_[i][2]) )
                idmap[tetIdx_[i][2]] = vtxCnt ++;
            if ( !idmap.count(tetIdx_[i][3]) )
                idmap[tetIdx_[i][3]] = vtxCnt ++;
        }
    }
}

template <typename T>
void TetMesh<T>::get_vtx_neighbor_table(std::vector< std::set<int> >& tbl) const
{
    tbl.resize(vertices_.size());
    for(size_t i = 0;i < vertices_.size();++ i) tbl[i].clear();

    for(size_t i = 0;i < tetIdx_.size();++ i)
    {
        tbl[tetIdx_[i][0]].insert(tetIdx_[i][1]);
        tbl[tetIdx_[i][0]].insert(tetIdx_[i][2]);
        tbl[tetIdx_[i][0]].insert(tetIdx_[i][3]);

        tbl[tetIdx_[i][1]].insert(tetIdx_[i][0]);
        tbl[tetIdx_[i][1]].insert(tetIdx_[i][2]);
        tbl[tetIdx_[i][1]].insert(tetIdx_[i][3]);

        tbl[tetIdx_[i][2]].insert(tetIdx_[i][0]);
        tbl[tetIdx_[i][2]].insert(tetIdx_[i][1]);
        tbl[tetIdx_[i][2]].insert(tetIdx_[i][3]);

        tbl[tetIdx_[i][3]].insert(tetIdx_[i][0]);
        tbl[tetIdx_[i][3]].insert(tetIdx_[i][1]);
        tbl[tetIdx_[i][3]].insert(tetIdx_[i][2]);
    }
}

#ifdef USE_NAMESPACE
}
#endif

#endif
