/******************************************************************************
 *  File: TetMeshWriter.hpp
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
#ifndef IO_TETMESH_WRITER_HPP
#   define IO_TETMESH_WRITER_HPP

#include "geometry/TetMesh.hpp"
#include "geometry/FixVtxTetMesh.hpp"
#include "utils/macros.h"
#include <fstream>

/*
 * write mesh into abaqus mesh format
 */
struct AbaqusMeshWriter
{
#ifdef USE_NAMESPACE
    template <typename T>
    static int write_mesh(const char* filename, const carbine::TetMesh<T>& mesh);
#else
    template <typename T>
    static int write_mesh(const char* filename, const TetMesh<T>& mesh);
#endif
};

#ifdef USE_NAMESPACE
template <typename T>
int AbaqusMeshWriter::write_mesh(const char* filename, const carbine::TetMesh<T>& mesh)
#else
template <typename T>
int AbaqusMeshWriter::write_mesh(const char* filename, const TetMesh<T>& mesh)
#endif
{
    using namespace std;
    ofstream fout(filename);
    if ( fout.fail() )
    {
        fprintf(stderr, "File %s cannot be open\n", filename);
        return ERROR_RETURN;
    }

    const std::vector< Point3<T> >& vtx = mesh.vertices();
    const std::vector< typename TetMesh<T>::TetIdx >& tdx = mesh.tet_indices();
    fout << "*NODE, NSET=NALL" << endl;
    fout << setprecision(20);
    for(size_t i = 0;i < vtx.size();++ i)
        fout << "    " << i+1 << ", "
             << vtx[i].x << ", "
             << vtx[i].y << ", "
             << vtx[i].z << endl;
    fout << "*ELEMENT, TYPE=C3D4, ELSET=EALL" << endl;
    for(size_t i = 0;i < tdx.size();++ i)
        fout << "    " << i+1 << ", "
             << tdx[i][0]+1 << ", " << tdx[i][1]+1 << ", "
             << tdx[i][2]+1 << ", " << tdx[i][3]+1 << endl;
    fout.close();
    return SUCC_RETURN;
}

struct FV_AbaqusMeshWriter
{
#ifdef USE_NAMESPACE
    template <typename T>
    static int write_mesh(const char* filename, const carbine::FixVtxTetMesh<T>& mesh);
#else
    template <typename T>
    static int write_mesh(const char* filename, const FixVtxTetMesh<T>& mesh);
#endif
};

#ifdef USE_NAMESPACE
template <typename T>
int FV_AbaqusMeshWriter::write_mesh(const char* filename, const carbine::FixVtxTetMesh<T>& mesh)
#else
template <typename T>
int FV_AbaqusMeshWriter::write_mesh(const char* filename, const FixVtxTetMesh<T>& mesh)
#endif
{
    using namespace std;
    ofstream fout(filename);
    if ( fout.fail() ) 
    {
        fprintf(stderr, "Cannot open file %s to write\n", filename);
        return ERROR_RETURN;
    }

    const vector< Point3<T> >& vtx = mesh.vertices();
    const vector< typename TetMesh<T>::TetIdx >& tdx = mesh.tet_indices();
    fout << "*NODE, NSET=NALL" << endl;
    fout << setprecision(20);
    for(size_t i = 0;i < vtx.size();++ i)
        fout << "    " << i+1 << ", "
             << vtx[i].x << ", "
             << vtx[i].y << ", "
             << vtx[i].z << endl;
    fout << "*ELEMENT, TYPE=C3D4, ELSET=EALL" << endl;
    for(size_t i = 0;i < tdx.size();++ i)
        fout << "    " << i+1 << ", "
             << tdx[i][0]+1 << ", " << tdx[i][1]+1 << ", "
             << tdx[i][2]+1 << ", " << tdx[i][3]+1 << endl;

    size_t nfixed = mesh.num_fixed_vertices();
    if ( nfixed > 0 ) 
    {
        fout << "*NSET, NSET=NFIXED" << endl;
        for(size_t i = 1;i <= nfixed;++ i)
        {
            if ( i % 10 == 1 ) fout << "    ";
            fout << i << ", ";
            if ( i % 10 == 0 ) fout << endl;
        }
        if ( nfixed % 10 != 0 ) fout << endl;
        fout << "*BOUNDARY" << endl
             << "    NFIXED,1,3"
             << endl;
    }

    fout.close();
    return SUCC_RETURN;
}

// ----------------------------------------------------------------------------

struct TetMeshWriter_Double
{
#ifdef USE_NAMESPACE
    template <typename T>
    static int write_mesh(const char* filename, const carbine::TetMesh<T>& mesh);
#else
    template <typename T>
    static int write_mesh(const char* filename, const TetMesh<T>& mesh);
#endif
};

#ifdef USE_NAMESPACE
template <typename T>
int TetMeshWriter_Double::write_mesh(const char* filename, const carbine::TetMesh<T>& mesh)
#else
template <typename T>
int TetMeshWriter_Double::write_mesh(const char* filename, const TetMesh<T>& mesh)
#endif
{
#ifdef USE_NAMESPACE
    using namespace carbine;
#endif

    FILE* file = fopen(filename, "wb");
    if ( !file )
    {
        fprintf(stderr, "File %s cannot be created\n", filename);
        return ERROR_RETURN;
    }

    size_t nb;

    const std::vector< Point3<T> >& vertices = mesh.rest_positions();
    int numVtx = (int)vertices.size();

    // output vertex array size
    nb = fwrite((void *)&numVtx, sizeof(int), 1, file);
    assert(nb == 1);
    // write the whole vertex array
    nb = fwrite((void *)&vertices[0], sizeof(Point3<T>), numVtx, file);
    if ( (int)nb != numVtx )
    {
        std::cerr << "ERROR: writing error to file: " << file << std::endl;
        exit(1);
    }

    const std::vector< typename TetMesh<T>::TetIdx >& indices = mesh.tet_indices();
    int numTet = (int)indices.size();
    // write the tet indices
    nb = fwrite((void *)&numTet, sizeof(int), 1, file);
    assert(nb == 1);
    
    for(int i = 0;i < numTet;++ i)
    {
        nb = fwrite((void *)&indices[i], sizeof(typename TetMesh<T>::TetIdx), 1, file);
        assert(nb == 1);
    }
    fclose(file);
    return SUCC_RETURN;
}

struct FV_TetMeshWriter_Double
{
#ifdef USE_NAMESPACE
    template <typename T>
    static int write_mesh(const char* filename, const carbine::FixVtxTetMesh<T>& mesh);

    template <typename T>
    static int write_mesh(const char* filename, const carbine::TetMesh<T>& mesh);
#else
    template <typename T>
    static int write_mesh(const char* filename, const FixVtxTetMesh<T>& mesh);

    template <typename T>
    static int write_mesh(const char* filename, const TetMesh<T>& mesh);
#endif
};

#ifdef USE_NAMESPACE
template <typename T>
int FV_TetMeshWriter_Double::write_mesh(const char* filename, 
        const carbine::TetMesh<T>& mesh)
#else
template <typename T>
int FV_TetMeshWriter_Double::write_mesh(const char* filename, 
        const TetMesh<T>& mesh)
#endif
{
#ifdef USE_NAMESPACE
    using namespace carbine;
#endif

    FILE* file = fopen(filename, "wb");
    if ( !file )
    {
        fprintf(stderr, "File %s cannot be created\n", filename);
        return ERROR_RETURN;
    }

    size_t nb;

    const std::vector< Point3<T> >& vertices = mesh.rest_positions();
    int numFree  = vertices.size();
    int numFixed = 0;

    // output vertex array size
    nb = fwrite((void *)&numFixed, sizeof(int), 1, file);
    assert(nb == 1);
    nb = fwrite((void *)&numFree, sizeof(int), 1, file);
    assert(nb == 1);

    // No fixed vertices, do nothing
    // write the free vertices
    nb = fwrite((void *)&vertices[numFixed], sizeof(Point3<T>), numFree, file);
    if ( (int)nb != numFree )
    {
        std::cerr << "ERROR: writing error to file: " << file << std::endl;
        exit(1);
    }

    const std::vector< typename TetMesh<T>::TetIdx >& indices = mesh.tet_indices();
    int numTet = (int)indices.size();
    // write the tet indices
    nb = fwrite((void *)&numTet, sizeof(int), 1, file);
    assert(nb == 1);
    
    for(int i = 0;i < numTet;++ i)
    {
        nb = fwrite((void *)&indices[i], sizeof(typename TetMesh<T>::TetIdx), 1, file);
        assert(nb == 1);
    }
    fclose(file);
    return SUCC_RETURN;
}

#ifdef USE_NAMESPACE
template <typename T>
int FV_TetMeshWriter_Double::write_mesh(const char* filename, 
        const carbine::FixVtxTetMesh<T>& mesh)
#else
template <typename T>
int FV_TetMeshWriter_Double::write_mesh(
        const char* filename, const FixVtxTetMesh<T>& mesh)
#endif
{
#ifdef USE_NAMESPACE
    using namespace carbine;
#endif

    FILE* file = fopen(filename, "wb");
    if ( !file )
    {
        fprintf(stderr, "File %s cannot be created\n", filename);
        return ERROR_RETURN;
    }

    size_t nb;

    const std::vector< Point3<T> >& vertices = mesh.rest_positions();
    int numFree  = mesh.num_free_vertices();
    int numFixed = mesh.num_fixed_vertices();

    // output vertex array size
    nb = fwrite((void *)&numFixed, sizeof(int), 1, file);
    assert(nb == 1);
    nb = fwrite((void *)&numFree, sizeof(int), 1, file);
    assert(nb == 1);

    // write the fixed vertices
    nb = fwrite((void *)&vertices[0], sizeof(Point3<T>), numFixed, file);
    if ( (int)nb != numFixed )
    {
        std::cerr << "ERROR: writing error to file: " << file << std::endl;
        exit(1);
    }
    // write the free vertices
    nb = fwrite((void *)&vertices[numFixed], sizeof(Point3<T>), numFree, file);
    if ( (int)nb != numFree )
    {
        std::cerr << "ERROR: writing error to file: " << file << std::endl;
        exit(1);
    }

    const std::vector< typename TetMesh<T>::TetIdx >& indices = mesh.tet_indices();
    int numTet = (int)indices.size();
    // write the tet indices
    nb = fwrite((void *)&numTet, sizeof(int), 1, file);
    assert(nb == 1);
    
    for(int i = 0;i < numTet;++ i)
    {
        nb = fwrite((void *)&indices[i], sizeof(typename TetMesh<T>::TetIdx), 1, file);
        assert(nb == 1);
    }
    fclose(file);
    return SUCC_RETURN;
}

#endif
