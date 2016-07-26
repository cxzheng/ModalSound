/******************************************************************************
 *  File: RigidObjRecorder.cpp
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
#include "RigidObjRecorder.h"
#include <fstream>
#include <iomanip>
#include <string>
#include "deformable/fem.hpp"
#include "io/MatrixIO.hpp"
#include "io/TglMeshWriter.hpp"
#include "utils/macros.h"

using namespace std;

namespace RigidObjRecorder
{

/*
 * save the geometry of the mesh in plain text
 * NOTE: all the index for the vertices into both surface mesh and tet mesh
 * are 0-based
 *
 * The format is
 * <# of surface vtx>
 * vid1_in_tet_mesh vid1_in_tgl_mesh n1.x n1.y n1.z a1  (for vertex 1)
 * vid2_in_tet_mesh vid1_in_tgl_mesh n2.x n2.y n2.z a2  (for vertex 2)
 * vid3_in_tet_mesh vid1_in_tgl_mesh n3.x n3.y n3.z a3  (for vertex 3)
 * ...
 */
int save_geometry(const TRigidBody* body, const char* file)
{
    ofstream fout(file);
    if ( fout.fail() ) 
    {
        cerr << "ERROR: Cannot open file: " << file << endl;
        return ERROR_RETURN;
    }

    map<int, int>  idmap;
    body->mesh()->surface_id_map(idmap);    // map from vertex id in tet mesh to that in surf mesh
    map<int, int>::iterator end = idmap.end();

    fout << idmap.size() << endl;
    fout << setprecision(22);
    const vector< Vector3<REAL> >& nml = body->collision_processor()->vtx_pseudo_nml();
    const valarray< REAL >&       area = body->collision_processor()->surface_mesh()->vertex_areas();
    for(map<int, int>::iterator it = idmap.begin();it != end;++ it)
    {
        Vector3<REAL> n = nml[it->second];
        n.normalize();
        fout << it->first << ' ' << it->second << ' ' 
             << n.x << ' ' << n.y << ' ' << n.z << ' '
             << area[it->second] << endl;
    }
    fout.close();

    return SUCC_RETURN;
}

int save_displacement(double ts, const TetMesh<double>* mesh, const char* file)
{
    ofstream fout(file, ios::binary);
    if ( fout.fail() )
    {
        cerr << "ERROR: Cannot open file: " << file << endl;
        return ERROR_RETURN;
    }

    fout.write((char *)&ts, sizeof(double));

    const vector<Point3d>& vtx = mesh->vertices();
    const vector<Point3d>& rst = mesh->rest_positions();
    int n = (int)vtx.size();
    Vector3d disp;
    fout.write((char *)&n, sizeof(int));
    for(int i = 0;i < n;++ i)
    {
        disp = vtx[i] - rst[i];
        fout.write((char *)&disp, sizeof(Vector3d));
    }

    fout.close();
    return SUCC_RETURN;
}

/*
 * save the surface mesh into an OBJ mesh file
 */
int save_surface_mesh(const TRigidBody* body, const char* file)
{
    return MeshObjWriter::write(
            *(body->collision_processor()->surface_mesh()), 
            file);
}

int save_mass_mat(const TRigidBody* body, const char* file)
{
    PardisoMatrix<REAL> M;
    DeformableFEM::mass_mat(body->mesh(), M);

/*
    //// write in matlab format (for testing)
#if defined(DEBUG) | defined(_DEBUG)
    string filename(file);
    filename += ".txt";
    ofstream fout(filename.c_str());
    fout << setprecision(22);
    PardisoMatrixIO::write_matlab_sparse(M, fout);
    fout.close();
#endif
*/
    return PardisoMatrixIO::write_csc_format(M, file);
}

int save_stiff_mat(const TRigidBody* body, const TMaterial* pmaterial, const char* file)
{
    PardisoMatrix<REAL> K;
    DeformableFEM::stiffness_mat(body->mesh(), pmaterial, K);

/*
//// write in matlab format (for testing)
#if defined(DEBUG) | defined(_DEBUG)
    string filename(file);
    filename += ".txt";
    ofstream fout(filename.c_str());
    fout << setprecision(22);
    PardisoMatrixIO::write_matlab_sparse(K, fout);
    fout.close();
#endif
*/
    return PardisoMatrixIO::write_csc_format(K, file);
}

} // end namespace

