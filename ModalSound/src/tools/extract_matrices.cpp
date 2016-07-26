/*
 * Read the given tetrahedron mesh, and use that to extract the matrices
 * such that the stiffness matrix and mass matrix
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "utils/term_msg.h"
#include "geometry/TriangleMesh.hpp"
#include "geometry/FixedVtxTetMesh.hpp"
#include "io/TetMeshReader.hpp"
#include "io/TglMeshWriter.hpp"
#include "io/MatrixIO.hpp"
#include "deformable/stvk.h"
#include "deformable/fem.hpp"

using namespace std;

typedef FixedVtxTetMesh<double>               TMesh;

static TMesh                mesh;
static TriangleMesh<double> tglmesh;
static double youngModulus;
static double poissonRatio;
static int    outType = 0;
static string tetfile = "";
static int    matOutType = 0;

static vector<Vector3d> tglNml;
static vector<Vector3d> vtxNml;

static void usage(char* cmd)
{
    printf("Usage: %s -f <tet_file> -v -y <Young's modulus value> -p <Poisson's ratio> [-m -k]\n", cmd);
    printf("    -f <tet_file>       Specify the file name tetrahedron mesh, no extension is needed.\n");
    printf("                        The extension of the tetrahedron mesh should be .tet\n");
    printf("    -y <value>          Young's modulus value\n");
    printf("    -p <value>          Poisson's ratio\n");
    printf("    -m                  Output mass matrix, named as <tet_file>.mass.csc\n");
    printf("    -k                  Output stiffness matrix, named as <tet_file>.stiff.csc\n");
    printf("    -g                  Output geometry file, named as <tet_file>.geo.txt\n");
    printf("    -s                  Output the surface triangle mesh\n");
    printf("    -d <value>          Sparse matrix format type.\n");
    printf("                        0 (default): csc binary format\n");
    printf("                        1          : MATLAB sparse binary format\n");
    printf("\n");
}

static void parse_cmd(int argc, char* argv[])
{
    int opt;
    while ( (opt = getopt(argc, argv, "hf:y:p:d:mkgs")) != -1 )
    {
        switch (opt)
        {
            case 'h':
                usage(argv[0]);
                exit(0);
            case 'f':
                tetfile = optarg;
                break;
            case 'y':
                youngModulus = atof(optarg);
                break;
            case 'p':
                poissonRatio = atof(optarg);
                break;
            case 'd':
                matOutType = atoi(optarg);
                break;
            case 'm':
                outType |= 1;
                break;
            case 'k':
                outType |= 2;
                break;
            case 'g':
                outType |= 4;
                break;
            case 's':
                outType |= 8;
                break;
        }
    }
}

static void compute_tgl_pseudo_normals()
{
    const vector<Point3d>&  vtx = tglmesh.vertices();
    const vector<Tuple3ui>& tgl = tglmesh.surface_indices();

    tglNml.resize(tgl.size());
    for(size_t i = 0;i < tgl.size();++ i)
    {
        tglNml[i] = Triangle<double>::normal(
                vtx[tgl[i][0]], vtx[tgl[i][1]], vtx[tgl[i][2]]);
        if ( tglNml[i].length_sqr() < 1E-24 )
        {
            fprintf(stderr, "ERROR: triangle has zero area: %.30lf\n",
                    tglNml[i].length_sqr());
            exit(1);
        }
        tglNml[i].normalize();
    }
}

static void compute_vtx_pseudo_normals()
{
    const vector<Point3d>&  vtx = tglmesh.vertices();
    const vector<Tuple3ui>& tgl = tglmesh.surface_indices();

    vtxNml.resize(vtx.size());
    memset(&vtxNml[0], 0, sizeof(Point3d)*vtxNml.size());

    for(size_t i = 0;i < tgl.size();++ i)
    {
        const Vector3<double>& nml = tglNml[i];

        vtxNml[tgl[i][0]] += nml * Triangle<double>::angle(
              vtx[tgl[i][2]], vtx[tgl[i][0]], vtx[tgl[i][1]]);
        vtxNml[tgl[i][1]] += nml * Triangle<double>::angle(
              vtx[tgl[i][0]], vtx[tgl[i][1]], vtx[tgl[i][2]]);
        vtxNml[tgl[i][2]] += nml * Triangle<double>::angle(
              vtx[tgl[i][1]], vtx[tgl[i][2]], vtx[tgl[i][0]]);
    }
}

static void load_mesh()
{
    char fname[128];
    sprintf(fname, "%s.tet", tetfile.c_str());
    printf("Load mesh [%s] ... ", fname);
    FV_TetMeshLoader_Double::load_mesh(fname, mesh);
    printf(" [OK]\n");

    printf("Initialize mesh ... ");
    mesh.init();
    mesh.update_surface();

    if ( outType & 4 || outType & 8 )
        mesh.extract_surface(&tglmesh);
    if ( outType & 4 )
    {
        tglmesh.update_vertex_areas();
        compute_tgl_pseudo_normals();
        compute_vtx_pseudo_normals();
    }
    printf(" [OK]\n");
}

static void output_mass_matrix()
{
    char fname[128];

    PardisoMatrix<double> M;
    DeformableFEM::mass_mat(&mesh, M);

    switch ( matOutType )
    {
        case 0:
            sprintf(fname, "%s.mass.csc", tetfile.c_str());
            PardisoMatrixIO::write_csc_format(M, fname);
            break;
        case 1:
            sprintf(fname, "%s.mass.spm", tetfile.c_str());
            PardisoMatrixIO::write_matlab_sparse_bin(M, fname);
            break;
        default:
            SHOULD_NEVER_HAPPEN(-1);
    }
    
    //ofstream fout(fname);
    //PardisoMatrixIO::write_matlab_sparse(M, fout);
    //fout.close();
}

static void output_stiff_matrix()
{
    char fname[128];

    PardisoMatrix<REAL> K;
    StVKMaterial material(
            (youngModulus*poissonRatio)/((1.+poissonRatio)*(1.-2.*poissonRatio)),
            youngModulus / (2.*(1.+poissonRatio)));
    DeformableFEM::stiffness_mat(&mesh, &material, K);

    switch ( matOutType )
    {
        case 0:
            sprintf(fname, "%s.stiff.csc", tetfile.c_str());
            PardisoMatrixIO::write_csc_format(K, fname);
            break;
        case 1:
            sprintf(fname, "%s.stiff.spm", tetfile.c_str());
            PardisoMatrixIO::write_matlab_sparse_bin(K, fname);
            break;
        default:
            SHOULD_NEVER_HAPPEN(-1);
    }
    
    //ofstream fout(fname);
    //PardisoMatrixIO::write_matlab_sparse(K, fout);
    //fout.close();
}

static void output_surface_mesh()
{
    char fname[128];
    sprintf(fname, "%s.obj", tetfile.c_str());
    MeshObjWriter::write(tglmesh, fname);
}

static void output_geometry()
{
    char fname[128];
    sprintf(fname, "%s.geo.txt", tetfile.c_str());
    ofstream fout(fname);
    if ( fout.fail() )
    {
        cerr << "ERROR: Cannot open file: " << fname << endl;
        return;
    }

    map<int, int> idmap;
    mesh.surface_id_map(idmap);
    const map<int, int>::iterator end = idmap.end();
    const valarray<double>& area = tglmesh.vertex_areas();

    fout << idmap.size() << endl;
    fout << setprecision(20);
    for(map<int,int>::iterator it = idmap.begin();it != end;++ it)
    {
        Vector3d n = vtxNml[it->second];
        n.normalize();
        fout << it->first << ' ' << it->second << ' '
             << n.x << ' ' << n.y << ' ' << n.z << ' '
             << area[it->second] << endl;
    }
    fout.close();
}

int main(int argc, char* argv[])
{
    parse_cmd(argc, argv);
    if ( youngModulus <= 1E-8 || poissonRatio < 0 )
    {
        PRINT_ERROR("Invalid Young's modulus or Poisson ratio is given\n");
        exit(1);
    }
    if ( tetfile.empty() ) 
    {
        PRINT_ERROR("No tetrahedron file is given\n");
        exit(1);
    }
    if ( !(outType & 15) )
    {
        PRINT_WARNING("No output type specified\n");
        exit(0);
    }

    load_mesh();
    if ( outType & 1 ) output_mass_matrix();
    if ( outType & 2 ) output_stiff_matrix();
    if ( outType & 4 ) output_geometry();
    if ( outType & 8 ) output_surface_mesh();

    return 0;
}
