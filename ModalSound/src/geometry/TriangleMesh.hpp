#ifndef GEOMETRY_TRIANGLE_MESH_HPP
#   define GEOMETRY_TRIANGLE_MESH_HPP

#include <assert.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <valarray>
#include <limits>
#include <vector>
#include <set>
#include <unordered_map>

#include "Triangle.hpp"
#include "Point3.hpp"
#include "sc/Vector3.hpp"

#ifdef USE_NAMESPACE
namespace sploosh
{
#endif

//! Mesh made by a group of triangles
template <typename T>
class TriangleMesh
{
    public:
        TriangleMesh() { }
        TriangleMesh(const TriangleMesh& m):totalArea_(m.totalArea_),
                vertices_(m.vertices_), normals_(m.normals_),
                triangles_(m.triangles_), vtxAreas_(m.vtxAreas_)
        { }

        struct NeighborRec
        {
            size_t   num;    // how many neighbor faces (3 at most)
            uint32_t id[3];
        };

        void clear()
        {
            vertices_.clear();
            normals_.clear();
            triangles_.clear();
        }

        //! Return the ID of the added vertices
        int add_vertex(const Point3<T>& vtx)
        {
            vertices_.push_back(vtx);
            return vertices_.size() - 1;
        }

        //! Return the ID of added vertices
        int add_vertex_normal(const Point3<T>& vtx, const Vector3<T>& nml)
        {
            vertices_.push_back(vtx);
            normals_.push_back(nml);
            return vertices_.size() - 1;
        }

        //! Return the current number of triangles
        int add_triangle(unsigned int v0, unsigned int v1, unsigned int v2)
        {
            assert(v0 != v1 && v0 != v2 && v1 != v2);
            assert(v0 < vertices_.size() &&
                   v1 < vertices_.size() &&
                   v2 < vertices_.size());
            triangles_.push_back(Tuple3ui(v0, v1, v2));
            return triangles_.size();
        }

        int num_vertices() const { return vertices_.size(); }
        int num_triangles() const { return triangles_.size(); }

        std::vector< Point3<T> >& vertices() 
        {  return vertices_; }

        std::vector<Tuple3ui>& triangles()
        {  return triangles_; }

        const std::vector< Point3<T> >& vertices() const
        {  return vertices_; }

        const std::vector<Tuple3ui>& surface_indices() const
        {  return triangles_; }

        const std::vector<Tuple3ui>& triangles() const
        {  return triangles_; }

        const std::vector< Vector3<T> >& normals() const
        {  return normals_; }

        const std::valarray<T>& vertex_areas() const
        {  return vtxAreas_; }

        const Tuple3ui& triangle_ids(uint32_t tid) const
        {
            assert(tid < triangles_.size());
            return triangles_[tid];
        }

        const Point3<T>& vertex(uint32_t vid) const
        {  
            assert(vid < vertices_.size());
            return vertices_[vid];
        }

        const Point3<T>* vertex_ptr(uint32_t vid) const
        {  
            assert(vid < vertices_.size());
            return &(vertices_[vid]);
        }

        const Vector3<T>& normal(uint32_t vid) const
        {
            assert(vid < vertices_.size() && normals_.size() == vertices_.size());
            return normals_[vid];
        }

        const Vector3<T>* normal_ptr(uint32_t vid) const
        {
            assert(vid < vertices_.size() && normals_.size() == vertices_.size());
            return &(normals_[vid]);
        }

        bool has_normals() const 
        {
            return !normals_.empty() && vertices_.size() == normals_.size();
        }

        bool empty() const 
        {  return vertices_.empty(); }

        T total_area() const 
        {  return totalArea_; }
        T triangle_area(uint32_t tglId) const
        {   
            return Triangle<T>::area(
                    vertices_[triangles_[tglId].x],
                    vertices_[triangles_[tglId].y],
                    vertices_[triangles_[tglId].z]);
        }

        void bounding_box(Point3<T>& low, Point3<T>& up) const;

        void generate_normals();
        void generate_pseudo_normals();
        //! Save the current mesh to a file
        int save_mesh_txt(const char* file) const;
        int load_mesh_txt(const char* file);
        void update_vertex_areas();

        void translate_x(T dx);
        void translate_y(T dy);
        void translate_z(T dz);

        //! Get the triangle face neighbors; for each triangle, return
        //  a list of its neighbor triangles
        void get_face_neighborship(std::vector<NeighborRec>&) const;
        //! Get the table of neighbors of all vertices
        void get_vtx_neighborship(std::vector< std::set<int> >&) const;
        //! Get the adjacent triangles of each vertex
        void get_vtx_tgls(std::vector< std::set<int> >&) const;

    private:
        T                           totalArea_;
        std::vector< Point3<T> >    vertices_;
        std::vector< Vector3<T> >   normals_;      // vertex normals
        std::vector<Tuple3ui>       triangles_;    // indices of triangle vertices
        std::valarray<T>            vtxAreas_;     // area of each triangles
};

// ---------------------------------------------------------------------------------
template <typename T>
void TriangleMesh<T>::translate_x(T dx)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 20000) shared(dx)
    for(int i = 0;i < (int)vertices_.size();++ i)
#else
    for(size_t i = 0;i < vertices_.size();++ i)
#endif
        vertices_[i].x += dx;
}

template <typename T>
void TriangleMesh<T>::translate_y(T dy)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 20000) shared(dy)
    for(int i = 0;i < (int)vertices_.size();++ i)
#else
    for(size_t i = 0;i < vertices_.size();++ i)
#endif
        vertices_[i].y += dy;
}

template <typename T>
void TriangleMesh<T>::translate_z(T dz)
{
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 20000) shared(dz)
    for(int i = 0;i < (int)vertices_.size();++ i)
#else
    for(size_t i = 0;i < vertices_.size();++ i)
#endif
        vertices_[i].z += dz;
}

template <typename T>
void TriangleMesh<T>::get_face_neighborship(std::vector<NeighborRec>& neighbors) const
{
    std::unordered_map< int, std::unordered_map<int, int> > hash;
    int v0, v1;

    neighbors.resize(triangles_.size());
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 20000) shared(neighbors)
    for(int i = 0;i < (int)triangles_.size();++ i)
#else
    for(size_t i = 0;i < triangles_.size();++ i)
#endif
        neighbors[i].num = 0;

    for(size_t i = 0;i < triangles_.size();++ i)
    {
        // edge 0
        v0 = triangles_[i].x;
        v1 = triangles_[i].y;
        if ( v0 > v1 ) std::swap(v0, v1);   // v0 < v1
        if ( hash.count(v0) && hash[v0].count(v1) )
        {
            int fid = hash[v0][v1];
            neighbors[i].id[neighbors[i].num ++] = fid;
            neighbors[fid].id[neighbors[fid].num ++] = i;
        }
        else
            hash[v0][v1] = i;

        // edge 1
        v0 = triangles_[i].y;
        v1 = triangles_[i].z;
        if ( v0 > v1 ) std::swap(v0, v1);
        if ( hash.count(v0) && hash[v0].count(v1) )
        {
            int fid = hash[v0][v1];
            neighbors[i].id[neighbors[i].num ++] = fid;
            neighbors[fid].id[neighbors[fid].num ++] = i;
        }
        else
            hash[v0][v1] = i;

        // edge 2
        v0 = triangles_[i].z;
        v1 = triangles_[i].x;
        if ( v0 > v1 ) std::swap(v0, v1);
        if ( hash.count(v0) && hash[v0].count(v1) )
        {
            int fid = hash[v0][v1];
            neighbors[i].id[neighbors[i].num ++] = fid;
            neighbors[fid].id[neighbors[fid].num ++] = i;
        }
        else
            hash[v0][v1] = i;
    }
}

template <typename T>
void TriangleMesh<T>::update_vertex_areas()
{
    const T s = (T)1 / (T)3;
    totalArea_ = 0;
    vtxAreas_.resize(vertices_.size(), (T)0);
    for(int i = 0;i < (int)triangles_.size();++ i)
    {
        T area = Triangle<T>::area(
                vertices_[triangles_[i].x],
                vertices_[triangles_[i].y],
                vertices_[triangles_[i].z]); 
        totalArea_ += area;
        area *= s;

        for(int j = 0;j < 3;++ j)
            vtxAreas_[triangles_[i][j]] += area;
    }
}

template <typename T>
void TriangleMesh<T>::generate_pseudo_normals()
{
    std::vector< Vector3<T> > tglnmls(triangles_.size());
    for(size_t i = 0;i < triangles_.size();++ i)
    {
        tglnmls[i] = Triangle<T>::normal(
                vertices_[triangles_[i].x],
                vertices_[triangles_[i].y],
                vertices_[triangles_[i].z]);
        if ( tglnmls[i].length_sqr() < 1E-18 )
        {
            fprintf(stderr, "ERROR: triangle has zero area: %.18g\n",
                    tglnmls[i].length_sqr());
            exit(1);
        }
        tglnmls[i].normalize();
    }

    normals_.resize(vertices_.size());
    memset(&normals_[0], 0, sizeof(Vector3<T>)*vertices_.size());
    for(size_t i = 0;i < triangles_.size();++ i)
    {
        const Vector3<T>& nml = tglnmls[i];
        normals_[triangles_[i].x] += nml * Triangle<T>::angle(
                vertices_[triangles_[i].z],
                vertices_[triangles_[i].x],
                vertices_[triangles_[i].y]);
        normals_[triangles_[i].y] += nml * Triangle<T>::angle(
                vertices_[triangles_[i].x],
                vertices_[triangles_[i].y],
                vertices_[triangles_[i].z]);
        normals_[triangles_[i].z] += nml * Triangle<T>::angle(
                vertices_[triangles_[i].y],
                vertices_[triangles_[i].z],
                vertices_[triangles_[i].x]);
    }

    for(size_t i = 0;i < vertices_.size();++ i)
        normals_[i].normalize();
}

template <typename T>
void TriangleMesh<T>::generate_normals()
{
    vtxAreas_.resize(vertices_.size(), (T)0);
    normals_.resize(vertices_.size());
    memset(&normals_[0], 0, sizeof(Vector3<T>)*vertices_.size());
    totalArea_ = 0;

    for(int i = 0;i < (int)triangles_.size();++ i)
    {
        Point3<T>& p0 = vertices_[triangles_[i].x];
        Point3<T>& p1 = vertices_[triangles_[i].y];
        Point3<T>& p2 = vertices_[triangles_[i].z];
        Vector3<T> nml = (p1 - p0).crossProduct(p2 - p0) * 0.5;  // weight = area
        T area = nml.length();
        totalArea_ += area;

        for(int j = 0;j < 3;++ j)
        {
            vtxAreas_[triangles_[i][j]] += area;
            normals_[triangles_[i][j]] += nml;
        }
    }

    const T alpha = (T)1 / (T)3;
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) schedule(dynamic, 10000)
    for(int i = 0;i < (int)vertices_.size();++ i)
#else
    for(int i = 0;i < (int)vertices_.size();++ i)
#endif
    {
        normals_[i].normalize();
        vtxAreas_[i] *= alpha;
    }
}

/*!
 * number of vertices [N/F] number of triangles
 * v0.x v0.y v0.z
 * v1.x v1.y v1.z
 * ......
 * triangle0.v0 triangle0.v1 triangle0.v2
 * triangle1.v0 triangle1.v1 triangle1.v2
 * triangle2.v0 triangle2.v1 triangle2.v2
 * ......
 */
template <typename T>
int TriangleMesh<T>::save_mesh_txt(const char* file) const
{
    using namespace std;

    ofstream fout(file);
    if ( !fout.good() ) return -1;

    bool nml_out;
    fout << vertices_.size();
    if ( !vertices_.empty() && vertices_.size() == normals_.size() )
    {
        fout << " N ";
        nml_out = true;
    }
    else
    {
        fout << " F ";
        nml_out = false;
    }
    fout << triangles_.size() << endl;

    fout << setprecision(16);
    for(size_t i = 0;i < vertices_.size();++ i)
        fout << vertices_[i].x 
             << ' ' << vertices_[i].y 
             << ' ' << vertices_[i].z << endl;
    // output normals
    if ( nml_out )
        for(size_t i = 0;i < normals_.size();++ i)
            fout << normals_[i].x << ' ' 
                 << normals_[i].y << ' ' 
                 << normals_[i].z << endl;

    for(size_t i = 0;i < triangles_.size();++ i)
        fout << triangles_[i].x << ' ' 
             << triangles_[i].y << ' ' 
             << triangles_[i].z << endl;

    fout.flush();
    fout.close();
    return 0;
}

template <typename T>
int TriangleMesh<T>::load_mesh_txt(const char* file)
{
    using namespace std;

    ifstream fin(file);
    if ( !fin.good() ) return -1;

    int nvtx, ntgl;
    char C;
    fin >> nvtx >> C >> ntgl;
    clear();
    T a, b, c;
    unsigned int v0, v1, v2;
    for(int i = 0;i < nvtx;++ i)
    {
        fin >> a >> b >> c;
        vertices_.push_back(Point3<T>(a, b, c));
    }
    if ( fin.fail() ) return -1;
    if ( C == 'N' )
    {
        for(int i = 0;i < nvtx;++ i)
        {
            fin >> a >> b >> c;
            normals_.push_back(Vector3<T>(a, b, c));
        }
        if ( fin.fail() ) return -1;
    }
    for(int i = 0;i < ntgl;++ i)
    {
        fin >> v0 >> v1 >> v2;
        triangles_.push_back(Tuple3ui(v0, v1, v2));
    }
    if ( fin.fail() ) return -1;

    fin.close();
    return 0;
}

template <typename T>
void TriangleMesh<T>::bounding_box(Point3<T>& low, Point3<T>& up) const
{
    const T inf = std::numeric_limits<T>::infinity();
    low.set(inf, inf, inf);
    up.set(-inf, -inf, -inf);

    const typename std::vector< Point3<T> >::const_iterator end = vertices_.end();
    for(typename std::vector< Point3<T> >::const_iterator it = vertices_.begin();
            it != end;++ it)
    {
        low.x = min(low.x, it->x);
        low.y = min(low.y, it->y);
        low.z = min(low.z, it->z);

        up.x = max(up.x, it->x);
        up.y = max(up.y, it->y);
        up.z = max(up.z, it->z);
    }
}

template <typename T>
void TriangleMesh<T>::get_vtx_neighborship(std::vector< std::set<int> >& tbl) const
{
    tbl.resize(vertices_.size());
    for(size_t i = 0;i < vertices_.size();++ i) tbl[i].clear();

    const std::vector<Tuple3ui>::const_iterator end = triangles_.end();
    for(std::vector<Tuple3ui>::const_iterator it = triangles_.begin();
            it != end;++ it)
    {
        tbl[it->x].insert(it->y); tbl[it->x].insert(it->z);
        tbl[it->y].insert(it->x); tbl[it->y].insert(it->z);
        tbl[it->z].insert(it->x); tbl[it->z].insert(it->y);
    }
}

/*
 * Return a list of triangles adjacent to each vertex
 */
template <typename T>
void TriangleMesh<T>::get_vtx_tgls(std::vector< std::set<int> >& vtx_ts) const
{
    vtx_ts.resize(vertices_.size());
    for(size_t i = 0;i < triangles_.size();++ i)
    {
        vtx_ts[ triangles_[i].x ].insert(i);
        vtx_ts[ triangles_[i].y ].insert(i);
        vtx_ts[ triangles_[i].z ].insert(i);
    }
}

#ifdef USE_NAMESPACE
}
#endif

#endif
