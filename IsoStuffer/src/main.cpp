/******************************************************************************
 *  File: main.cpp
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
/*
 * This program is to generate the tetrahedron mesh from an isosurface.
 *
 * It basically implement the paper: isosurface stuffing: fast tetrahedral meshes
 * with good dihedral angles, siggraph 2007
 */
#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include "IsoStuffer.hpp"
#include "DistValProc.h"
#include "io/TglMeshReader.hpp"
#include "io/TetMeshWriter.hpp"

using namespace std;

typedef DistFunc<double, DistValProc>           TDistFunc;
typedef IsoStuffing::IsoStuffer<TDistFunc>      TIsoStuffer;

static int res     = 6;
static int nl      = 4;
static int nmargin = 7;
static double alphaL = 0.25;
static double alphaS = 0.42978;

static TriangleMesh<double>*    ptglmesh  = NULL;
static DistValProc*             dvproc    = NULL;
static TDistFunc*               pdistfunc = NULL;
static TIsoStuffer*             pstuffer  = NULL;
static TetMesh<double>*         pmesh     = NULL;

static void clean_up()
{
    delete pmesh;
    delete pstuffer;
    delete pdistfunc;
    delete dvproc;
    delete ptglmesh;
}

int main(int argc, char* argv[])
{
    namespace po = boost::program_options;
    po::options_description genericOpt("Generic options");
    genericOpt.add_options()
            ("help,h", "display help information");
    po::options_description configOpt("Configuration");
    configOpt.add_options()
        ("res,R",    po::value<int>(&res), "bounding box resolution")
        ("nlevel,L", po::value<int>(&nl),  "number of levels")
        ("margin,M", po::value<int>(&nmargin), "margin of the bounding box")
        ("alpha,a",  po::value<double>(&alphaL), "alpha long value in isostuff algorithm")
        ("beta,b",  po::value<double>(&alphaS), "alpha short value in isostuff algorithm");
    // use configure file to specify the option
    po::options_description cfileOpt("Configure file");
    cfileOpt.add_options()
        ("files", po::value< vector<string> >(), "input files")
        ("cfg-file", po::value<string>(), "configuration file");

    po::options_description cmdOpts;
    cmdOpts.add(genericOpt).add(configOpt).add(cfileOpt);

    po::positional_options_description p;
    p.add("files", 2);

    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdOpts).positional(p).run(), vm);
    if ( vm.count("cfg-file") )
    {
        ifstream ifs(vm["cfg-file"].as<string>().c_str());
        store(parse_config_file(ifs, configOpt), vm);
    }
    po::notify(vm);

    if ( vm.count("help") )
    {
        printf("Usage: %s [options] <input obj file> <output tet file>\n", argv[0]);
        cout << cmdOpts << endl;
        return 0;
    }
    if ( !vm.count("files") || vm["files"].as< vector<string> >().size() != 2 ) 
    {
        cerr << "No input/output files are given" << endl;
        printf("Usage: %s [options] <input obj file> <output tet file>\n", argv[0]);
        cout << cmdOpts << endl;
        return 1;
    }
    res = 1 << res;
    const vector<string>& fs = vm["files"].as< vector<string> >();
    cout << "==============================================" << endl;
    cout << "  res         = " << res << endl;
    cout << "  num level   = " << nl  << endl;
    cout << "  margin      = " << nmargin << endl;
    cout << "  alpha_Long  = " << alphaL << endl;
    cout << "  alpha_Short = " << alphaS << endl;
    cout << "  input/output: " << fs[0] << ' ' << fs[1] << endl;
    cout << "==============================================" << endl;

    ptglmesh = new TriangleMesh<double>;
    MeshObjReader::read(fs[0].c_str(), *ptglmesh);
    dvproc   = new DistValProc(ptglmesh, res, nl);
    pdistfunc= new TDistFunc(dvproc);

    dvproc->update_mesh_params(res, nl, nmargin);
    pstuffer = new TIsoStuffer(dvproc->resolution(), nl, dvproc->grid_size(),
            dvproc->min_point(), alphaL, alphaS, *pdistfunc);
    pmesh = new TetMesh<double>;
    pstuffer->create_mesh(pmesh);
    TetMeshWriter_Double::write_mesh(fs[1].c_str(), *pmesh);

    clean_up();
    return 0;
}
