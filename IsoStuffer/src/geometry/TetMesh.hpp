/******************************************************************************
 *  File: TetMesh.hpp
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
#ifndef GEOMETRY_TETMESH_HPP
#   define GEOMETRY_TETMESH_HPP

#include "config.h"
#include <vector>
#include <valarray>
#include <map>
#include <set>

#ifdef USE_HASH_MAP
#   include <unordered_map>
#else
#   include <tr1/unordered_map>
#endif

#include <string.h>
#include "Tet.hpp"
#include "TriangleMesh.hpp"
#include "utils/tuple.hpp"
#include "linearalgebra/Tuple4.hpp"

#ifdef USE_NAMESPACE
namespace carbine
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
        TetMesh():m_canAddVtx(true) { }
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
            memcpy(&m_vertices[0], &m_restPos[0], sizeof(Point3<T>)*m_restPos.size());
        }

        template <typename FromT>
        void add_vertex(const Point3<FromT>& v)
        {
#if defined(DEBUG) | defined(_DEBUG)
            assert(m_canAddVtx);
#endif

            m_vertices.push_back(v);
            m_restPos.push_back(v);
        }

        /*
         * add a vertex with returning the id of the newly added vertex
         */
        template <typename FromT>
        int add_vertex_rt(const Point3<FromT>& v)
        {
#if defined(DEBUG) | defined(_DEBUG)
            assert(m_canAddVtx);
#endif

            m_vertices.push_back(v);
            m_restPos.push_back(v);

            return m_vertices.size() - 1;
        }

        /*
         * ### NOTE: ###
         * Be aware of what you are doing when calling this method. It may 
         * screw up the pointers in m_tets
         *
         * return the ID of the added vertex
         */
        template <typename FromT>
        int add_vertex_unsafe(const Point3<FromT>& now, const Point3<FromT>& rest)
        {
            m_restPos.push_back(rest);
            m_vertices.push_back(now);
            return m_vertices.size() - 1;
        }

        void update_deformation_gradient();

        /*
         * If any new vertex is added after some tet has been added. The pointer in
         * m_vertices may point to wrong address. The method should be called to make
         * sure the pointer in Tet<T> maintains the correct address.
         *
         * void conform_vtx_pointers();
         */

        /* ============ retrival methods ============ */

        /*! return number of tetrahedrons */
        size_t num_tets() const { return m_tetIdx.size(); }
        size_t num_vertices() const { return m_vertices.size(); }
        size_t num_surface_tgls() const { return m_surfIdx.size(); }

        const std::vector<TetIdx>& tet_indices() const { return m_tetIdx; }

        const std::vector<Tuple3ui>& surface_indices() const
        { return m_surfIdx; }

        /*! return masses associated with each vertices */
        const std::vector<T>& masses() const 
        { return m_masses; }

        ///*! return inverse of masses associated with each vertices */
        //const std::vector<T>& inverse_masses() const
        //{ return m_invMasses; }

        const std::vector< Tet<T> >& tets() const 
        { return m_tets; }

        std::vector< Tet<T> >& tets() 
        { return m_tets; }

        std::vector< Point3<T> >& vertices()
        { return m_vertices; }

        const Point3<T>& vertex(size_t vid) const
        { 
            assert(vid < m_vertices.size());
            return m_vertices[vid]; 
        }

        const Point3<T>& rest_vertex(size_t vid) const
        {
            assert(vid < m_restPos.size());
            return m_restPos[vid]; 
        }

        const std::vector< Vector3<T> >& normals() const
        { return m_normals; }

        const std::vector< Point3<T> >& vertices() const
        { return m_vertices; }

        const std::vector< Point3<T> >& rest_positions() const
        { return m_restPos; }

        std::vector< Point3<T> >& rest_positions()
        { return m_restPos; }

#ifdef WITH_DEFORMABLE_TET
        std::vector< Vector3<T> >& velocities()
        { return m_vel; }

        std::vector< Vector3<T> >& accelerations()
        { return m_acc; }
#endif

        /*
         * m_surf2tet is computed when update_surface is called
         */
        const std::vector<int>& surf_to_tet() const
        { return m_surf2tet; }

        /* 
         * check if the mesh is write referenced
         * only for debug purpose
         */
        bool check_mesh() const
        {
            for(size_t i = 0;i < m_tetIdx.size();++ i)
            {
                for(int j = 0;j < 4;++ j)
                    if ( &(m_vertices[m_tetIdx[i][j]]) != m_tets[i].vtx(j) ) 
                    {
                        fprintf(stderr, "WARNING: CHECK Mesh Failed at tet [%d] %d\n", 
                                (int)i, (int)m_tetIdx[i][j]);
                        return false;
                    }
            }
            return true;
        }

    private:
        void add_normal(unsigned int a, unsigned int b, unsigned int c,
                std::vector<int>& cnts);

    protected:
        std::vector< Point3<T> >    m_vertices;     // current positions of vertices
        std::vector< Point3<T> >    m_restPos;      // position at rest pose

    private:
#if defined(DEBUG) | defined(_DEBUG)
        bool        m_canAddVtx;
#endif

#ifdef WITH_DEFORMABLE_TET
        std::vector< Vector3<T> >   m_vel;          // velocity of each vertex
        std::vector< Vector3<T> >   m_acc;          // acceleration of each vertex
#endif
        std::vector< Tet<T> >       m_tets;
        std::vector< TetIdx >       m_tetIdx;       // vertex index for each vertex of the tet
        std::vector< Tuple3ui >     m_surfIdx;      // vertices' index for surface triangles
        /*!
         * the vector<int> which maps each triangle in m_surfIdx to
         * tet that has the triangle
         */
        std::vector<int>            m_surf2tet;

        std::vector<T>              m_masses;       // m_masses[i] is the mass of m_vertices[i]
        //std::vector<T>              m_invMasses;    // m_invMasses[i] is the inverse of mass of m_vertices[i]
        std::vector< Vector3<T> >   m_normals;
};

///////////////////////////////////////////////////////////////////////////////

template <typename T>
TetMesh<T>::TetMesh(const TetMesh<T>& rhs)
{
    clear();
    for(size_t i = 0;i < rhs.m_vertices.size();++ i)
        add_vertex_unsafe(rhs.m_vertices[i], rhs.m_restPos[i]);

    for(size_t i = 0;i < rhs.m_tets.size();++ i)
        add_tet(rhs.m_tetIdx[i][0], rhs.m_tetIdx[i][1], 
                rhs.m_tetIdx[i][2], rhs.m_tetIdx[i][3]);

#ifdef WITH_DEFORMABLE_TET
    //// initialize it
    m_vel.resize(m_vertices.size());
    m_acc.resize(m_vertices.size());
    memset(&m_vel[0], 0, sizeof(Vector3<T>)*m_vel.size());
    memset(&m_acc[0], 0, sizeof(Vector3<T>)*m_acc.size());
#endif

    m_masses.resize(m_vertices.size());
    //m_invMasses.resize(m_vertices.size());
    m_tets.resize(m_tetIdx.size());
    for(size_t i = 0;i < m_tetIdx.size();++ i)
    {
        m_tets[i].init(&m_restPos[m_tetIdx[i][0]], &m_restPos[m_tetIdx[i][1]],
                       &m_restPos[m_tetIdx[i][2]], &m_restPos[m_tetIdx[i][3]]);
        T v = m_tets[i].volume() * 0.25;

        m_tets[i].conform_vtx_pointers(
                &m_vertices[m_tetIdx[i][0]],
                &m_vertices[m_tetIdx[i][1]],
                &m_vertices[m_tetIdx[i][2]],
                &m_vertices[m_tetIdx[i][3]]);

        m_masses[m_tetIdx[i][0]] += v;
        m_masses[m_tetIdx[i][1]] += v;
        m_masses[m_tetIdx[i][2]] += v;
        m_masses[m_tetIdx[i][3]] += v;
    }
#if defined(DEBUG) | defined(_DEBUG)
    m_canAddVtx = false;
#endif
}

template <typename T>
void TetMesh<T>::uninitialized_clone(const TetMesh<T>* rhs)
{
    clear();

    m_vertices.resize(rhs->m_vertices.size());
    m_restPos.resize(rhs->m_restPos.size());
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(rhs)
#endif
    for(int i = 0;i < (int)rhs->m_vertices.size();++ i)
    {
        m_vertices[i] = rhs->m_vertices[i];
        m_restPos[i]  = rhs->m_restPos[i];
    }

    m_tetIdx.resize(rhs->m_tetIdx.size());
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(rhs)
#endif
    for(int i = 0;i < (int)rhs->m_tetIdx.size();++ i)
        m_tetIdx[i] = rhs->m_tetIdx[i];

#ifdef WITH_DEFORMABLE_TET
    //// initialize it
    m_vel.resize(m_vertices.size());
    m_acc.resize(m_vertices.size());
    memset(&m_vel[0], 0, sizeof(Vector3<T>)*m_vel.size());
    memset(&m_acc[0], 0, sizeof(Vector3<T>)*m_acc.size());
#endif

    m_surfIdx.resize(rhs->m_surfIdx.size());
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(rhs)
#endif
    for(int i = 0;i < (int)rhs->m_surfIdx.size();++ i)
        m_surfIdx[i] = rhs->m_surfIdx[i];

    if ( !rhs->m_surf2tet.empty() )
    {
        m_surf2tet.resize(rhs->m_surf2tet.size());
#ifdef USE_OPENMP
        #pragma omp parallel for default(none) shared(rhs)
#endif
        for(int i = 0;i < (int)rhs->m_surf2tet.size();++ i)
            m_surf2tet[i] = rhs->m_surf2tet[i];
    }

    if ( !rhs->m_normals.empty() )
    {
        m_normals.resize(rhs->m_normals.size());
#ifdef USE_OPENMP
        #pragma omp parallel for default(none) shared(rhs)
#endif
        for(int i = 0;i < (int)rhs->m_normals.size();++ i)
            m_normals[i] = rhs->m_normals[i];
    }
}

/*!
 * - initialize the velocity and acceleration vector
 * - compute the mass and inverse mass vector
 *
 * NOTE: we create the m_tets array here, assuming that no add_vertex is gonna call
 *       after this. Otherwise, a growing m_vertex array may change its address, and 
 *       hence screw up the vertex pointers in m_tet
 */
template <typename T>
void TetMesh<T>::init()
{
#ifdef WITH_DEFORMABLE_TET
    m_vel.resize(m_vertices.size());
    m_acc.resize(m_vertices.size());
    memset(&m_vel[0], 0, sizeof(Vector3<T>)*m_vel.size());
    memset(&m_acc[0], 0, sizeof(Vector3<T>)*m_acc.size());
#endif

    m_masses.resize(m_vertices.size());
    //m_invMasses.resize(m_vertices.size());
    m_tets.resize(m_tetIdx.size());
    for(size_t i = 0;i < m_tetIdx.size();++ i)
    {
        m_tets[i].init(&m_vertices[m_tetIdx[i][0]],
                &m_vertices[m_tetIdx[i][1]],
                &m_vertices[m_tetIdx[i][2]],
                &m_vertices[m_tetIdx[i][3]]);
        T v = m_tets[i].volume() * 0.25;

        m_masses[m_tetIdx[i][0]] += v;
        m_masses[m_tetIdx[i][1]] += v;
        m_masses[m_tetIdx[i][2]] += v;
        m_masses[m_tetIdx[i][3]] += v;
    }
#if defined(DEBUG) | defined(_DEBUG)
    m_canAddVtx = false;
#endif

    /*
    for(size_t i = 0;i < m_masses.size();++ i)
        if ( m_masses[i] < 1E-20 )
        {
            fprintf(stderr, "ERROR: error on Tetrahedron %d with ZERO volumn %.24lf\n", 
                    (int)i, m_masses[i]);
            m_invMasses[i] = 0;
        }
        else
            m_invMasses[i] = 1. / m_masses[i];
    */
}

/*
 * NOTE: don't call 
 *  m_tets.push_back(Tet<T>(&m_vertices[v0], &m_vertices[v1],
 *              &m_vertices[v2], &m_vertices[v3]));
 * here. This is because m_vertices is still growing. System might 
 * change the address of m_vertices
 */
template <typename T>
void TetMesh<T>::add_tet(int v0, int v1, int v2, int v3)
{
    assert(v0 < (int)m_vertices.size() && v1 < (int)m_vertices.size() &&
           v2 < (int)m_vertices.size() && v3 < (int)m_vertices.size());

    //m_tets.push_back(Tet<T>(&m_vertices[v0], &m_vertices[v1],
    //            &m_vertices[v2], &m_vertices[v3]));
    m_tetIdx.push_back(TetIdx(v0, v1, v2, v3));
}

template <typename T>
void TetMesh<T>::clear()
{
    m_vertices.clear();
    m_restPos.clear();
#ifdef WITH_DEFORMABLE_TET    
    m_vel.clear();
    m_acc.clear();
#endif
    m_tets.clear();
    m_tetIdx.clear();
    m_surfIdx.clear();
    m_masses.clear();
    //m_invMasses.clear();
    m_normals.clear();
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
#ifdef USE_UNORDERED_MAP
    using namespace std::tr1;
#else
    using namespace std;
#endif
    unordered_map< int, unordered_map<int, unordered_map<int, int> > > hash;

    int v0, v1, v2;

    neighbors.resize(m_tetIdx.size());
    for(size_t i = 0;i < m_tetIdx.size(); ++ i)
        neighbors[i].num = 0;

    for(size_t i = 0;i < m_tetIdx.size();++ i)
    {
        if ( i < tetRemoved.size() && tetRemoved[i] ) continue;

        // face 0
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][3];
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
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][2];
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
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
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
        v0 = m_tetIdx[i][1];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
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
#ifdef USE_UNORDERED_MAP
    using namespace std::tr1;
#else
    using namespace std;
#endif
    unordered_map< int, unordered_map<int, unordered_map<int, int> > > hash;

    int v0, v1, v2;

    neighbors.resize(m_tetIdx.size());
    for(size_t i = 0;i < m_tetIdx.size(); ++ i)
        neighbors[i].num = 0;

    for(size_t i = 0;i < m_tetIdx.size();++ i)
    {
        // face 0
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][3];
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
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][2];
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
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
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
        v0 = m_tetIdx[i][1];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
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
#ifdef USE_UNORDERED_MAP
    using namespace std::tr1;
#else
    using namespace std;
#endif
    unordered_map< int, unordered_map< int, unordered_map<int, int> > > hash;
    unordered_map< int, unordered_map< int, unordered_map<int, int> > > hash2;

    int v0, v1, v2;

    if ( m_normals.size() != m_vertices.size() )
        m_normals.resize(m_vertices.size());

#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 2000)
    for(size_t i = 0;i < m_vertices.size();++ i)
#else
    for(int i = 0;i < (int)m_vertices.size();++ i)
#endif
        m_normals[i].zero();

    for(size_t i = 0;i < m_tetIdx.size();++ i)
    {
        // face 0
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
        hash2[v0][v1][v2] = i;

        // face 1
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][2];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
        hash2[v0][v1][v2] = i;

        // face 2
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
        hash2[v0][v1][v2] = i;

        // face 3
        v0 = m_tetIdx[i][1];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
        hash2[v0][v1][v2] = i;
    }

    m_surfIdx.clear();
    m_surf2tet.clear();

    std::vector<int> cnts(m_vertices.size());
    memset(&cnts[0], 0, sizeof(int)*cnts.size());

    for(size_t i = 0;i < m_tetIdx.size();++ i)
    {
        // face 0
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            m_surfIdx.push_back(Tuple3i(
                        m_tetIdx[i][0], 
                        m_tetIdx[i][1],
                        m_tetIdx[i][3]));
            m_surf2tet.push_back(hash2[v0][v1][v2]);
            add_normal(m_tetIdx[i][0], m_tetIdx[i][1], m_tetIdx[i][3], cnts);
        }

        // face 1
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][2];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            m_surfIdx.push_back(Tuple3i(
                        m_tetIdx[i][0], 
                        m_tetIdx[i][2],
                        m_tetIdx[i][1]));
            m_surf2tet.push_back(hash2[v0][v1][v2]);
            add_normal(m_tetIdx[i][0], m_tetIdx[i][2], m_tetIdx[i][1], cnts);
        }

        // face 2
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            m_surfIdx.push_back(Tuple3i(
                        m_tetIdx[i][3], 
                        m_tetIdx[i][2],
                        m_tetIdx[i][0]));
            m_surf2tet.push_back(hash2[v0][v1][v2]);
            add_normal(m_tetIdx[i][3], m_tetIdx[i][2], m_tetIdx[i][0], cnts);
        }

        // face 3
        v0 = m_tetIdx[i][1];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            m_surfIdx.push_back(Tuple3i(
                        m_tetIdx[i][1], 
                        m_tetIdx[i][2],
                        m_tetIdx[i][3]));
            m_surf2tet.push_back(hash2[v0][v1][v2]);
            add_normal(m_tetIdx[i][1], m_tetIdx[i][2], m_tetIdx[i][3], cnts);
        }
    }

#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 2000) shared(cnts)
    for(size_t i = 0;i < m_vertices.size();++ i)
#else
    for(int i = 0;i < (int)m_vertices.size();++ i)
#endif
    {
        m_normals[i] /= (T)cnts[i];
        m_normals[i].normalize();
    }
    printf("Update surface: %d triangles\n", (int)m_surfIdx.size());
}

template <typename T>
void TetMesh<T>::update_deformation_gradient()
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 2000)
    for(int i = 0;i < (int)m_tets.size();++ i)
#else
    for(size_t i = 0;i < m_tets.size();++ i)
#endif
        m_tets[i].deformation_gradient();
}

template <typename T>
void TetMesh<T>::add_normal(unsigned int a, unsigned int b, 
        unsigned int c, std::vector<int>& cnts)
{
    Vector3<T> n = (m_restPos[b] - m_restPos[a]).crossProduct(
                m_restPos[c] - m_restPos[a]);

    ++ cnts[a];
    ++ cnts[b];
    ++ cnts[c];
    m_normals[a] += n;
    m_normals[b] += n;
    m_normals[c] += n;
}

template <typename T>
void TetMesh<T>::extract_surface(TriangleMesh<T>* smesh) const
{
#ifdef USE_UNORDERED_MAP
    using namespace std::tr1;
#else
    using namespace std;
#endif
    unordered_map< int, unordered_map< int, unordered_map<int, int> > > hash;
    unordered_map< int, int > idmap;

    smesh->clear();

    int v0, v1, v2;
    for(size_t i = 0;i < m_tetIdx.size();++ i)
    {
        // face 0
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 1
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][2];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 2
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 3
        v0 = m_tetIdx[i][1];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
    }

    for(size_t i = 0;i < m_tetIdx.size();++ i)
    {
        // face 0
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(m_tetIdx[i][0]) )
                idmap[m_tetIdx[i][0]] = smesh->add_vertex(m_restPos[m_tetIdx[i][0]]);
            if ( !idmap.count(m_tetIdx[i][1]) )
                idmap[m_tetIdx[i][1]] = smesh->add_vertex(m_restPos[m_tetIdx[i][1]]);
            if ( !idmap.count(m_tetIdx[i][3]) )
                idmap[m_tetIdx[i][3]] = smesh->add_vertex(m_restPos[m_tetIdx[i][3]]);
            smesh->add_triangle(idmap[m_tetIdx[i][0]], idmap[m_tetIdx[i][1]], idmap[m_tetIdx[i][3]]);
        }

        // face 1
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][2];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(m_tetIdx[i][0]) )
                idmap[m_tetIdx[i][0]] = smesh->add_vertex(m_restPos[m_tetIdx[i][0]]);
            if ( !idmap.count(m_tetIdx[i][2]) )
                idmap[m_tetIdx[i][2]] = smesh->add_vertex(m_restPos[m_tetIdx[i][2]]);
            if ( !idmap.count(m_tetIdx[i][1]) )
                idmap[m_tetIdx[i][1]] = smesh->add_vertex(m_restPos[m_tetIdx[i][1]]);
            smesh->add_triangle(idmap[m_tetIdx[i][0]], idmap[m_tetIdx[i][2]], idmap[m_tetIdx[i][1]]);
        }

        // face 2
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(m_tetIdx[i][0]) )
                idmap[m_tetIdx[i][0]] = smesh->add_vertex(m_restPos[m_tetIdx[i][0]]);
            if ( !idmap.count(m_tetIdx[i][3]) )
                idmap[m_tetIdx[i][3]] = smesh->add_vertex(m_restPos[m_tetIdx[i][3]]);
            if ( !idmap.count(m_tetIdx[i][2]) )
                idmap[m_tetIdx[i][2]] = smesh->add_vertex(m_restPos[m_tetIdx[i][2]]);
            smesh->add_triangle(idmap[m_tetIdx[i][0]], idmap[m_tetIdx[i][3]], idmap[m_tetIdx[i][2]]);
        }

        // face 3
        v0 = m_tetIdx[i][1];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(m_tetIdx[i][1]) )
                idmap[m_tetIdx[i][1]] = smesh->add_vertex(m_restPos[m_tetIdx[i][1]]);
            if ( !idmap.count(m_tetIdx[i][2]) )
                idmap[m_tetIdx[i][2]] = smesh->add_vertex(m_restPos[m_tetIdx[i][2]]);
            if ( !idmap.count(m_tetIdx[i][3]) )
                idmap[m_tetIdx[i][3]] = smesh->add_vertex(m_restPos[m_tetIdx[i][3]]);
            smesh->add_triangle(idmap[m_tetIdx[i][1]], idmap[m_tetIdx[i][2]], idmap[m_tetIdx[i][3]]);
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
#ifdef USE_UNORDERED_MAP
    using namespace std::tr1;
#else
    using namespace std;
#endif
    unordered_map< int, unordered_map< int, unordered_map<int, int> > > hash;

    idmap.clear();
    int v0, v1, v2;
    for(size_t i = 0;i < m_tetIdx.size();++ i)
    {
        // face 0
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 1
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][2];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 2
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];

        // face 3
        v0 = m_tetIdx[i][1];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        ++ hash[v0][v1][v2];
    }

    int vtxCnt = 0;
    for(size_t i = 0;i < m_tetIdx.size();++ i)
    {
        // face 0
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(m_tetIdx[i][0]) )
                idmap[m_tetIdx[i][0]] = vtxCnt ++;
            if ( !idmap.count(m_tetIdx[i][1]) )
                idmap[m_tetIdx[i][1]] = vtxCnt ++;
            if ( !idmap.count(m_tetIdx[i][3]) )
                idmap[m_tetIdx[i][3]] = vtxCnt ++;
        }

        // face 1
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][1];
        v2 = m_tetIdx[i][2];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(m_tetIdx[i][0]) )
                idmap[m_tetIdx[i][0]] = vtxCnt ++;
            if ( !idmap.count(m_tetIdx[i][2]) )
                idmap[m_tetIdx[i][2]] = vtxCnt ++;
            if ( !idmap.count(m_tetIdx[i][1]) )
                idmap[m_tetIdx[i][1]] = vtxCnt ++;
        }

        // face 2
        v0 = m_tetIdx[i][0];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(m_tetIdx[i][0]) )
                idmap[m_tetIdx[i][0]] = vtxCnt ++;
            if ( !idmap.count(m_tetIdx[i][3]) )
                idmap[m_tetIdx[i][3]] = vtxCnt ++;
            if ( !idmap.count(m_tetIdx[i][2]) )
                idmap[m_tetIdx[i][2]] = vtxCnt ++;
        }

        // face 3
        v0 = m_tetIdx[i][1];
        v1 = m_tetIdx[i][2];
        v2 = m_tetIdx[i][3];
        sort_triple(v0, v1, v2);
        if ( hash[v0][v1][v2] == 1 )
        {
            if ( !idmap.count(m_tetIdx[i][1]) )
                idmap[m_tetIdx[i][1]] = vtxCnt ++;
            if ( !idmap.count(m_tetIdx[i][2]) )
                idmap[m_tetIdx[i][2]] = vtxCnt ++;
            if ( !idmap.count(m_tetIdx[i][3]) )
                idmap[m_tetIdx[i][3]] = vtxCnt ++;
        }
    }
}

/*
 * construct a table: tbl[i] maintains the id of tets associated with the 
 * i-th vertex
 */
template <typename T>
void TetMesh<T>::get_vtx_neighbor_table(std::vector< std::set<int> >& tbl) const
{
    tbl.resize(m_vertices.size());
    for(size_t i = 0;i < m_vertices.size();++ i) tbl[i].clear();

    for(size_t i = 0;i < m_tetIdx.size();++ i)
    {
        tbl[m_tetIdx[i][0]].insert(m_tetIdx[i][1]);
        tbl[m_tetIdx[i][0]].insert(m_tetIdx[i][2]);
        tbl[m_tetIdx[i][0]].insert(m_tetIdx[i][3]);

        tbl[m_tetIdx[i][1]].insert(m_tetIdx[i][0]);
        tbl[m_tetIdx[i][1]].insert(m_tetIdx[i][2]);
        tbl[m_tetIdx[i][1]].insert(m_tetIdx[i][3]);

        tbl[m_tetIdx[i][2]].insert(m_tetIdx[i][0]);
        tbl[m_tetIdx[i][2]].insert(m_tetIdx[i][1]);
        tbl[m_tetIdx[i][2]].insert(m_tetIdx[i][3]);

        tbl[m_tetIdx[i][3]].insert(m_tetIdx[i][0]);
        tbl[m_tetIdx[i][3]].insert(m_tetIdx[i][1]);
        tbl[m_tetIdx[i][3]].insert(m_tetIdx[i][2]);
    }
}

#ifdef USE_NAMESPACE
}
#endif

#endif
