/******************************************************************************
 *  File: TglMeshReader.hpp
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
#ifndef GEOMETRIC_MESH_READER
#   define GEOMETRIC_MESH_READER

#include <fstream>
#include <vector>
#include "geometry/TriangleMesh.hpp"
#include "utils/macros.h"
#include "utils/strings.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

class MeshObjReader
{
    public:
        static const std::string  SUB_DELIM;

        template <typename T>
        static int read(const char* file, TriangleMesh<T>& mesh, 
                bool centerize = false, bool reversetglrot = false)
        {
            using namespace std;

            ifstream fin(file);
            if ( fin.fail() ) return ERROR_RETURN;
            
            char text[1024];
            vector<Point3d>     vtx;
            vector<Vector3d>    nml;
            vector<Tuple3ui>    tgl;

            int l = 1;
            fin.getline(text, 1024);
            while ( !fin.fail() )
            {
                vector<string> tokens = carbine::tokenize(string(text));
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
                        cerr << "(A) Incorrect file format at Line " << l << "!" << endl;
                        exit(1);
                    }
                    vtx.push_back(Point3d(carbine::Double(tokens[1]),
                                          carbine::Double(tokens[2]),
                                          carbine::Double(tokens[3])));
                }
                else if ( tokens[0] == "vn" )
                {
                    if ( tokens.size() != 4 )
                    {
                        cerr << "(B) Incorrect file format at Line " << l << "!" << endl;
                        exit(1);
                    }
                    nml.push_back(Vector3d(carbine::Double(tokens[1]),
                                           carbine::Double(tokens[2]),
                                           carbine::Double(tokens[3])));
                    nml.back().normalize();
                }
                else if ( tokens[0] == "f" )
                {
                    if ( tokens.size() != 4 )
                    {
                        cerr << "(C) Incorrect file format at Line " << l << "!" << endl;
                        exit(1);
                    }

                    vector<string> ts0 = carbine::tokenize(tokens[1], SUB_DELIM);
                    vector<string> ts1 = carbine::tokenize(tokens[2], SUB_DELIM);
                    vector<string> ts2 = carbine::tokenize(tokens[3], SUB_DELIM);

                    if ( ts0.size() != ts1.size() || ts1.size() != ts2.size() )
                    {
                        cerr << "(D) Incorrect file format at Line " << l << "!" << endl;
                        exit(1);
                    }

                    if ( ts0.size() > 1 )
                    {
                        if ( ts0[0] != ts0[1] || ts1[0] != ts1[1] || ts2[0] != ts2[1] )
                            cerr << "Warinning: Inconsistent file content at Line " << l << "!" << endl;

                        tgl.push_back(Tuple3ui(carbine::Int(ts0[0])-1,
                                               carbine::Int(ts1[0])-1,
                                               carbine::Int(ts2[0])-1));
                    }
                    else
                    {
                        //// the index in the OBJ file is 1-based
                        tgl.push_back(Tuple3ui(carbine::Int(tokens[1])-1,
                                               carbine::Int(tokens[2])-1,
                                               carbine::Int(tokens[3])-1));
                    }
                }

                fin.getline(text, 1024);
                ++ l;
            }

            fin.close();

            if ( !nml.empty() && nml.size() != vtx.size() )
            {
                cerr << "The file should include the same number of normals and vertices!" << endl;
                exit(1);
            }

            /* No triangles at all */
            if ( tgl.empty() )
                cerr << "Warning: THERE IS NO TRIANGLE MESHS AT ALL!" << endl;

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
                for(size_t i = 0;i < nml.size();++ i)
                    mesh.add_vertex_normal(vtx[i], nml[i]);
            }

            if ( reversetglrot )
                for(size_t i = 0;i < tgl.size();++ i)
                    mesh.add_triangle(tgl[i][0], tgl[i][2], tgl[i][1]);
            else
                for(size_t i = 0;i < tgl.size();++ i)
                    mesh.add_triangle(tgl[i][0], tgl[i][1], tgl[i][2]);

            return SUCC_RETURN;
        }

        /*!
         * Read the obj file for which the normals, if specified, should be indicated
         * at each vertex.
         */
        template <typename T>
        static int read_with_per_vtx_nml(const char* file, TriangleMesh<T>& mesh,
                                         bool reversetglrot = false)
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
                vector<string> tokens = carbine::tokenize(string(text));
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
                        cerr << "(a) Incorrect file format at Line " << l << "!" << endl;
                        exit(1);
                    }
                    vtx.push_back(Point3d(carbine::Double(tokens[1]),
                                          carbine::Double(tokens[2]),
                                          carbine::Double(tokens[3])));
                }
                else if ( tokens[0] == "vn" )
                {
                    if ( tokens.size() != 4 )
                    {
                        cerr << "(b) Incorrect file format at Line " << l << "!" << endl;
                        exit(1);
                    }
                    nml.push_back(Vector3d(carbine::Double(tokens[1]),
                                           carbine::Double(tokens[2]),
                                           carbine::Double(tokens[3])));
                    nml.back().normalize();
                }
                else if ( tokens[0] == "f" )
                {
                    if ( tokens.size() != 4 )
                    {
                        cerr << "(c) Incorrect file format at Line " << l << "!" << endl;
                        exit(1);
                    }

                    vector<string> ts0 = carbine::split(tokens[1], SUB_DELIM[0]);
                    vector<string> ts1 = carbine::split(tokens[2], SUB_DELIM[0]);
                    vector<string> ts2 = carbine::split(tokens[3], SUB_DELIM[0]);

                    if ( ts0.size() != ts1.size() || ts1.size() != ts2.size() )
                    {
                        cerr << "(d) Incorrect file format at Line " << l << "!" << endl;
                        exit(1);
                    }

                    tgl.push_back(Tuple3ui(carbine::Int(ts0[0])-1,
                                           carbine::Int(ts1[0])-1,
                                           carbine::Int(ts2[0])-1));
                    if ( ts0.size() > 2 )
                        tNml.push_back(Tuple3ui(carbine::Int(ts0[2])-1,
                                                carbine::Int(ts1[2])-1,
                                                carbine::Int(ts2[2])-1));
                }

                fin.getline(text, 1024);
                ++ l;
            }

            fin.close();

            if ( !nml.empty() && nml.size() != vtx.size() )
            {
                cerr << "The file should include the same number of normals and vertices!" << endl;
                exit(1);
            }

            if ( tgl.empty() ) cerr << "Warning: THERE IS NO TRIANGLE MESHS AT ALL!" << endl;

            if ( nml.empty() )
            {
                for(size_t i = 0;i < vtx.size();++ i)
                    mesh.add_vertex(vtx[i]);
            }
            else
            {
                vector<int> nmlmap(vtx.size(), -1);
                for(size_t i = 0;i < tgl.size();++ i)
                {
                    if ( nmlmap[tgl[i][0]] >= 0 && nmlmap[tgl[i][0]] != tNml[i][0] )
                    {
                        cerr << "inconsistent vertex normal map" << endl;
                        exit(1);
                    }
                    else
                        nmlmap[tgl[i][0]] = tNml[i][0];

                    if ( nmlmap[tgl[i][1]] >= 0 && nmlmap[tgl[i][1]] != tNml[i][1] )
                    {
                        cerr << "inconsistent vertex normal map" << endl;
                        exit(1);
                    }
                    else
                        nmlmap[tgl[i][1]] = tNml[i][1];

                    if ( nmlmap[tgl[i][2]] >= 0 && nmlmap[tgl[i][2]] != tNml[i][2] )
                    {
                        cerr << "inconsistent vertex normal map" << endl;
                        exit(1);
                    }
                    else
                        nmlmap[tgl[i][2]] = tNml[i][2];
                }

                for(size_t i = 0;i < nml.size();++ i)
                    mesh.add_vertex_normal(vtx[i], nml[nmlmap[i]]);
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
