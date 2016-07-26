/******************************************************************************
 *  File: StellarIO.hpp
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
#ifndef STELLAR_IO_HPP
#   define STELLAR_IO_HPP

#include <fstream>
#include "geometry/FixVtxTetMesh.hpp"
#include "utils/macros.h"

struct StellarTetMeshLoader
{
#ifdef USE_NAMESPACE
    template <typename T>
    static int load_mesh(const char* filename, carbine::TetMesh<T>* mesh)
#else
    template <typename T>
    static int load_mesh(const char* filename, TetMesh<T>* mesh)
#endif
    {
        using namespace std;
        char fname[128];

        sprintf(fname, "%s.node", filename);
        ifstream fin(fname);
        if ( fin.fail() )
        {
            fprintf(stderr, "Cannot open file %s to read\n", fname);
            return ERROR_RETURN;
        }

        int npt, nf, nattr, mark;
        fin >> npt >> nf >> nattr >> mark;
        if ( nf != 3 )
        {
            fprintf(stderr, "Incorrect format (NF != 3)\n");
            return ERROR_RETURN;
        }

        T xx, yy, zz;
        int initid, dummy;
        int totF = nattr + mark;
        // read vertex positions
        for(int rep = 0;rep < npt;++ rep)
        {
            if  ( !rep ) 
                fin >> initid;
            else 
                fin >> dummy;
            fin >> xx >> yy >> zz;
            for(int j = 0;j < totF;++ j) fin >> dummy;
            mesh->add_vertex(Point3d(xx, yy, zz));
        }
        fin.close();

        // read tetrahedrons
        sprintf(fname, "%s.ele", filename);
        fin.open(fname);
        if ( fin.fail() ) 
        {
            fprintf(stderr, "Cannot open file %s to read\n", fname);
            return ERROR_RETURN;
        }
        fin >> npt >> nf >> nattr;
        if ( nf != 4 || nattr )
        {
            fprintf(stderr, "Incorrect format (NF != 4)\n");
            return ERROR_RETURN;
        }
        int vid[4];
        for(int rep = 0;rep < npt;++ rep)
        {
            fin >> dummy >> vid[0] >> vid[1] >> vid[3] >> vid[2];
            mesh->add_tet(vid[0]-initid, vid[1]-initid, vid[2]-initid, vid[3]-initid);
        }
        fin.close();

        return SUCC_RETURN;
    }
};

#endif
