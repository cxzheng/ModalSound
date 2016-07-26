#ifndef IO_TETMESH_LOADER_HPP
#   define IO_TETMESH_LOADER_HPP

#include "geometry/FixedVtxTetMesh.hpp"
#include "utils/macros.h"

///////////////////////////////////////////////////////////////////////////////
//// TetMeshLoader_Double

/*! 
 * load tetrahedra mesh, all the data are in double precision (double) 
 *
 * It assumes there is no constraint vertices in the tet tetrahedron mesh.
 *
 * Use FV_TetMeshLoader_Double below to read the tet mesh with some fixed
 * vertices
 */
struct TetMeshLoader_Double
{
#ifdef USE_NAMESPACE
    template <typename T>
    static int load_mesh(const char* filename, sploosh::TetMesh<T>& mesh);
    template <typename T>
    static int load_mesh(const char* filename, sploosh::TetMesh<T>* mesh);
    template <typename T>
    static int load_mesh(const char* filename, sploosh::TetMesh<T>* mesh, 
                         T scale, const Vector3<T>& move);
#else
    template <typename T>
    static int load_mesh(const char* filename, TetMesh<T>& mesh);
    template <typename T>
    static int load_mesh(const char* filename, TetMesh<T>* mesh);
    template <typename T>
    static int load_mesh(const char* filename, TetMesh<T>* mesh,
                         T scale, const Vector3<T>& move);
#endif
};

/*!
 * NOTE: we don't initialize the mesh (we don't call mesh.init()) in this method.
 */
#ifdef USE_NAMESPACE
template <typename T>
int TetMeshLoader_Double::load_mesh(const char* filename, sploosh::TetMesh<T>& mesh)
#else
template <typename T>
int TetMeshLoader_Double::load_mesh(const char* filename, TetMesh<T>& mesh)
#endif
{
    FILE* file = fopen(filename, "rb");
    if ( !file )
    {
        fprintf(stderr, "File %s not found\n", filename);
        return ERROR_RETURN;
    }

    mesh.clear();

    int numVtx;
#ifndef NDEBUG
    size_t nb;
    // output vertex array sizes
    nb = fread((void*)&(numVtx), sizeof(int), 1, file);
    assert(nb == 1);
#else
    fread((void*)&(numVtx), sizeof(int), 1, file);
#endif

    printf("Read tetrahedron mesh:\n");
    printf(" %i vertices\n", numVtx);

    Point3d node;
    // read vertex positions
    for(int x = 0; x < numVtx; ++ x)
    {
#ifndef NDEBUG
        nb = fread((void*)&node, sizeof(double), 3, file);
        assert(nb == 3);
#else
        fread((void*)&node, sizeof(double), 3, file);
#endif
        mesh.add_vertex(node);
    }

    // output tet vertex lists
    int totalTets;
    int idx[4];
#ifndef NDEBUG
    nb = fread((void*)&totalTets, sizeof(int), 1, file);
    assert(nb == 1);
#else
    fread((void*)&totalTets, sizeof(int), 1, file);
#endif

    printf(" %i total tets\n", totalTets);
    for(int x = 0; x < totalTets;++ x)
    {
#ifndef NDEBUG
        nb = fread((void*)idx, sizeof(int), 4, file);
        assert(nb == 4);
#else
        fread((void*)idx, sizeof(int), 4, file);
#endif
        mesh.add_tet(idx[0], idx[1], idx[2], idx[3]);
    }

    fclose(file);
    return SUCC_RETURN;
}

/*!
 * NOTE: we don't initialize the mesh (we don't call mesh.init()) in this method.
 */
#ifdef USE_NAMESPACE
template <typename T>
int TetMeshLoader_Double::load_mesh(const char* filename, sploosh::TetMesh<T>* mesh)
#else
template <typename T>
int TetMeshLoader_Double::load_mesh(const char* filename, TetMesh<T>* mesh)
#endif
{
    FILE* file = fopen(filename, "rb");
    if ( !file )
    {
        fprintf(stderr, "File %s not found\n", filename);
        return ERROR_RETURN;
    }

    mesh->clear();

    int numVtx;
#ifndef NDEBUG
    size_t nb;
    // output vertex array sizes
    nb = fread((void*)&(numVtx), sizeof(int), 1, file);
    assert(nb == 1);
#else
    fread((void*)&(numVtx), sizeof(int), 1, file);
#endif

    printf("Read tetrahedron mesh:\n");
    printf(" %i vertices\n", numVtx);

    Point3d node;
    // read vertex positions
    for(int x = 0; x < numVtx; ++ x)
    {
#ifndef NDEBUG
        nb = fread((void*)&node, sizeof(double), 3, file);
        assert(nb == 3);
#else
        fread((void*)&node, sizeof(double), 3, file);
#endif
        mesh->add_vertex(node);
    }

    // output tet vertex lists
    int totalTets;
    int idx[4];
#ifndef NDEBUG
    nb = fread((void*)&totalTets, sizeof(int), 1, file);
    assert(nb == 1);
#else
    fread((void*)&totalTets, sizeof(int), 1, file);
#endif

    printf(" %i total tets\n", totalTets);
    for(int x = 0; x < totalTets;++ x)
    {
#ifndef NDEBUG
        nb = fread((void*)idx, sizeof(int), 4, file);
        assert(nb == 4);
#else
        fread((void*)idx, sizeof(int), 4, file);
#endif
        mesh->add_tet(idx[0], idx[1], idx[2], idx[3]);
    }

    fclose(file);
    return SUCC_RETURN;
}

#ifdef USE_NAMESPACE
template <typename T>
int TetMeshLoader_Double::load_mesh(const char* filename, 
        sploosh::TetMesh<T>* mesh, T scale, 
        const Vector3<T>& move)
#else
template <typename T>
int TetMeshLoader_Double::load_mesh(const char* filename, 
        TetMesh<T>* mesh, T scale, const Vector3<T>& move)
#endif
{
    FILE* file = fopen(filename, "rb");
    if ( !file )
    {
        fprintf(stderr, "File %s not found\n", filename);
        return ERROR_RETURN;
    }

    mesh->clear();

    int numVtx;
#ifndef NDEBUG
    size_t nb;
    // output vertex array sizes
    nb = fread((void*)&(numVtx), sizeof(int), 1, file);
    assert(nb == 1);
#else
    fread((void*)&(numVtx), sizeof(int), 1, file);
#endif

    printf("Read tetrahedron mesh:\n");
    printf(" %i vertices\n", numVtx);

    Point3d minPt(1E+10, 1E+10, 1E+10), maxPt(-1E+10, -1E+10, -1E+10);
    Point3d node;
    // read vertex positions
    for(int x = 0; x < numVtx; ++ x)
    {
#ifndef NDEBUG
        nb = fread((void*)&node, sizeof(double), 3, file);
        assert(nb == 3);
#else
        fread((void*)&node, sizeof(double), 3, file);
#endif
        node *= scale;
        node += move;
        maxPt.x = fmax(maxPt.x, node.x);
        maxPt.y = fmax(maxPt.y, node.y);
        maxPt.z = fmax(maxPt.z, node.z);
        minPt.x = fmin(minPt.x, node.x);
        minPt.y = fmin(minPt.y, node.y);
        minPt.z = fmin(minPt.z, node.z);
        mesh->add_vertex(node);
    }
    printf("  The object bounding box: [%.8lf %.8lf %lf] -> [%.8lf %.8lf %.8lf]\n",
            minPt.x, minPt.y, minPt.z, 
            maxPt.x, maxPt.y, maxPt.z);

    // output tet vertex lists
    int totalTets;
    int idx[4];
#ifndef NDEBUG
    nb = fread((void*)&totalTets, sizeof(int), 1, file);
    assert(nb == 1);
#else
    fread((void*)&totalTets, sizeof(int), 1, file);
#endif

    printf(" %i total tets\n", totalTets);
    for(int x = 0; x < totalTets;++ x)
    {
#ifndef NDEBUG
        nb = fread((void*)idx, sizeof(int), 4, file);
        assert(nb == 4);
#else
        fread((void*)idx, sizeof(int), 4, file);
#endif
        mesh->add_tet(idx[0], idx[1], idx[2], idx[3]);
    }

    fclose(file);
    return SUCC_RETURN;
}

///////////////////////////////////////////////////////////////////////////////
/*
 * Read mesh with Fixed vertices
 */
struct FV_TetMeshLoader_Double
{
#ifdef USE_NAMESPACE
    template <typename T>
    static int load_mesh(const char* filename, sploosh::FixedVtxTetMesh<T>& mesh);
#else
    template <typename T>
    static int load_mesh(const char* filename, FixedVtxTetMesh<T>& mesh);
#endif
};

#ifdef USE_NAMESPACE
template <typename T>
int FV_TetMeshLoader_Double<T>::load_mesh(
        const char* filename, sploosh::FixedVtxTetMesh<T>& mesh)
#else
template <typename T>
int FV_TetMeshLoader_Double::load_mesh(
        const char* filename, FixedVtxTetMesh<T>& mesh)
#endif
{
    FILE* file = fopen(filename, "rb");
    if ( !file )
    {
        fprintf(stderr, "File %s not found\n", filename);
        return ERROR_RETURN;
    }

    mesh.clear();

    int freeNumVtx, fixedNumVtx;
#ifndef NDEBUG
    size_t nb;
    nb = fread((void*)&fixedNumVtx, sizeof(int), 1, file);
    assert(nb == 1);
    nb = fread((void*)&freeNumVtx, sizeof(int), 1, file);
    assert(nb == 1);
#else
    fread((void*)&fixedNumVtx, sizeof(int), 1, file);
    fread((void*)&freeNumVtx, sizeof(int), 1, file);
#endif

    printf("INFO: Read tetrahedron mesh: %d free vertices, %d fixed vertices\n",
            freeNumVtx, fixedNumVtx);
    Point3d vtx;
    // read fixed vertices
    for(int i = 0;i < fixedNumVtx;++ i)
    {
#ifndef NDEBUG
        nb = fread((void*)&vtx, sizeof(double), 3, file);
        assert(nb == 3);
#else
        fread((void*)&vtx, sizeof(double), 3, file);
#endif
        mesh.add_vertex(vtx);
    }

    // read unconstrained vertices
    for(int i = 0;i < freeNumVtx;++ i)
    {
#ifndef NDEBUG
        nb = fread((void*)&vtx, sizeof(double), 3, file);
        assert(nb == 3);
#else
        fread((void*)&vtx, sizeof(double), 3, file);
#endif
        mesh.add_vertex(vtx);
    }
    mesh.set_fixed_vtx(fixedNumVtx);

    // add tetrahedron
    int ntet;
    int idx[4];
#ifndef NDEBUG
    nb = fread((void *)&ntet, sizeof(int), 1, file);
    assert(nb == 1);
#else
    fread((void *)&ntet, sizeof(int), 1, file);
#endif
    printf("INFO: Read %d tets\n", ntet);
    for(int i = 0;i < ntet;++ i)
    {
#ifndef NDEBUG
        nb = fread((void *)idx, sizeof(int), 4, file);
        assert(nb == 4);
#else
        fread((void *)idx, sizeof(int), 4, file);
#endif
        mesh.add_tet(idx[0], idx[1], idx[2], idx[3]);
    }
    fclose(file);
    return SUCC_RETURN;
}

#endif

