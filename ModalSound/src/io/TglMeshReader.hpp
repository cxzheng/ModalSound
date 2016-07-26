#ifndef GEOMETRIC_MESH_READER
#   define GEOMETRIC_MESH_READER

#include <fstream>
#include <vector>
#include "geometry/TriangleMesh.hpp"
#include "utils/macros.h"
#include "utils/term_msg.h"
#include "utils/strings.hpp"

#ifdef USE_NAMESPACE
namespace sploosh
{
#endif

class MeshObjReader
{
    public:
        /*!
         * Read the obj file for which the normals, if specified, should be indicated
         * at each vertex. In other words, the same vertex at different faces should
         * be specified with the same normal.
         */
        template <typename T>
        static int read(const char* file, TriangleMesh<T>& mesh, 
                bool centerize = false, bool reversetglrot = false)
        {
            using namespace std;

            ifstream fin(file);
            if ( fin.fail() ) return ERROR_RETURN;
            
            char text[1024];
            vector<Point3d>     vtx;    // temporary space to put the vertex
            vector<Vector3d>    nml;
            vector<Tuple3ui>    tgl;
            vector<Tuple3ui>    tNml;

            int l = 1;
            fin.getline(text, 1024);
            while ( !fin.fail() )
            {
                vector<string> tokens = sploosh::tokenize(string(text));
                if ( tokens.empty() ) 
                {
                    fin.getline(text, 1024);
                    ++ l;
                    continue;
                }

                if ( tokens[0] == "v" )
                {
                    if ( tokens.size() != 4 )
                    {
                        PRINT_ERROR("(A) Incorrect file format at Line %d.\n", l);
                        exit(1);
                    }
                    vtx.push_back(Point3d(sploosh::Double(tokens[1]),
                                          sploosh::Double(tokens[2]),
                                          sploosh::Double(tokens[3])));
                }
                else if ( tokens[0] == "vn" )
                {
                    if ( tokens.size() != 4 )
                    {
                        PRINT_ERROR("(B) Incorrect file format at Line %d.\n", l);
                        exit(1);
                    }
                    nml.push_back(Vector3d(sploosh::Double(tokens[1]),
                                           sploosh::Double(tokens[2]),
                                           sploosh::Double(tokens[3])));
                    nml.back().normalize();
                }
                else if ( tokens[0] == "f" )
                {
                    if ( tokens.size() != 4 )
                    {
                        PRINT_ERROR("(C) Incorrect file format at Line %d (# of fields is larger than 3)\n", l);
                        exit(1);
                    }

                    vector<string> ts0 = sploosh::split(tokens[1], '/');
                    vector<string> ts1 = sploosh::split(tokens[2], '/');
                    vector<string> ts2 = sploosh::split(tokens[3], '/');

                    if ( ts0.size() != ts1.size() || ts1.size() != ts2.size() )
                    {
                        PRINT_ERROR("(D) Incorrect file format at Line %d.\n", l);
                        exit(1);
                    }

                    tgl.push_back(Tuple3ui(sploosh::Int(ts0[0])-1,
                                           sploosh::Int(ts1[0])-1,
                                           sploosh::Int(ts2[0])-1));
                    if ( ts0.size() > 2 )
                        tNml.push_back(Tuple3ui(sploosh::Int(ts0[2])-1,
                                                sploosh::Int(ts1[2])-1,
                                                sploosh::Int(ts2[2])-1));
                }

                fin.getline(text, 1024);
                ++ l;
            }

            fin.close();

            if ( !nml.empty() && nml.size() != vtx.size() )
            {
                PRINT_ERROR("The file should include the same number of normals and vertices\n");
                exit(1);
            }

            /* No triangles at all */
            if ( tgl.empty() ) PRINT_WARNING("THERE IS NO TRIANGLE MESHS AT ALL!\n");

            /* Centerize the objects */
            if ( centerize )
            {
                Tuple3d maxPt(-1E+100, -1E+100, -1E+100), 
                        minPt(1E+100, 1E+100, 1E+100);
                for(size_t i = 0;i < vtx.size();++ i)
                {
                    maxPt.x = fmax(vtx[i].x, maxPt.x);
                    maxPt.y = fmax(vtx[i].y, maxPt.y);
                    maxPt.z = fmax(vtx[i].z, maxPt.z);

                    minPt.x = fmin(vtx[i].x, minPt.x);
                    minPt.y = fmin(vtx[i].y, minPt.y);
                    minPt.z = fmin(vtx[i].z, minPt.z);
                }

                Point3d center = (maxPt + minPt) * 0.5;
                for(size_t i = 0;i < vtx.size();++ i)
                    vtx[i] -= center;
            }

            if ( nml.empty() )
            {
                for(size_t i = 0;i < vtx.size();++ i)
                    mesh.add_vertex(vtx[i]);
            }
            else
            {
                int* nmlmap = new int[vtx.size()];
                memset(nmlmap, 0xFF, sizeof(int)*vtx.size());

                for(size_t i = 0;i < tgl.size();++ i)
                {
                    if ( nmlmap[tgl[i][0]] >= 0 && nmlmap[tgl[i][0]] != tNml[i][0] )
                    {
                        PRINT_ERROR("inconsistent vertex normal reference at triangle #%d.\n", (int)i);
                        delete[] nmlmap;
                        exit(1);
                    }
                    else
                        nmlmap[tgl[i][0]] = tNml[i][0];

                    if ( nmlmap[tgl[i][1]] >= 0 && nmlmap[tgl[i][1]] != tNml[i][1] )
                    {
                        PRINT_ERROR("inconsistent vertex normal reference at triangle #%d.\n", (int)i);
                        delete[] nmlmap;
                        exit(1);
                    }
                    else
                        nmlmap[tgl[i][1]] = tNml[i][1];

                    if ( nmlmap[tgl[i][2]] >= 0 && nmlmap[tgl[i][2]] != tNml[i][2] )
                    {
                        PRINT_ERROR("inconsistent vertex normal reference at triangle #%d.\n", (int)i);
                        delete[] nmlmap;
                        exit(1);
                    }
                    else
                        nmlmap[tgl[i][2]] = tNml[i][2];
                }

                for(size_t i = 0;i < nml.size();++ i)
                    mesh.add_vertex_normal(vtx[i], nml[nmlmap[i]]);
                delete[] nmlmap;
            }

            if ( reversetglrot )
                for(size_t i = 0;i < tgl.size();++ i)
                    mesh.add_triangle(tgl[i][0], tgl[i][2], tgl[i][1]);
            else
                for(size_t i = 0;i < tgl.size();++ i)
                    mesh.add_triangle(tgl[i][0], tgl[i][1], tgl[i][2]);

            return SUCC_RETURN;
        }
};

#ifdef USE_NAMESPACE
}
#endif
#endif
