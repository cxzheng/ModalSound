/******************************************************************************
 *  File: ellip_gen.cpp
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
 * generate the tetrahedron of the ellipsoid
 *
 * The levelset function is 
 * f(x,y,z) = x^2/a^2 + y^2/b^2 + z^2/c^2 - 1
 * 
 */
#include <iostream>

#include "IsoStuffer.hpp"
#include "io/TetMeshWriter.hpp"

using namespace std;

class EllipsoidDistFunc
{
    private:
        double invA2;   // 1 / a^2
        double invB2;   // 1 / b^2
        double invC2;   // 1 / c^2

    public:
        EllipsoidDistFunc(double a, double b, double c):
                invA2(1./(a*a)), invB2(1./(b*b)), invC2(1./(c*c))
        { }

        inline double operator() (const Point3d& pt) const
        {
            //return pt.norm() - 1.00123456789;
            return pt.x*pt.x*invA2 + pt.y*pt.y*invB2 + pt.z*pt.z*invC2 - 1.0;
        }
};

typedef IsoStuffing::IsoStuffer<EllipsoidDistFunc>  TIsoStuffer;

static const int MARGIN = 10;
static const int RES    = 1 << 6;
static const int NLEVEL = 6;
static const double ALPHA_LONG  = 0.25; 
static const double ALPHA_SHORT = 0.42978;

static TIsoStuffer*     pstuffer = NULL;
static TetMesh<double>* pmesh = NULL;

static void clean_up()
{
    delete pstuffer;
    delete pmesh;
}

int main(int argc, char* argv[])
{
    if ( argc != 5 )
    {
        cerr << "Usage: " << argv[0] << " a b c <out_mesh>" << endl;
        cerr << "       a b c: are the principle components of inertia tensor" << endl;
        exit(1);
    }

    const double Ix = atof(argv[1]); 
    const double Iy = atof(argv[2]); 
    const double Iz = atof(argv[3]);
    if ( Iy + Iz - Ix < 0 || Iz + Ix - Iy < 0 ||
         Ix + Iy - Iz < 0 )
    {
        cerr << "ERROR on given inertia values" << endl;
        return 2;
    }
    const double hd[3] = {sqrt((Iy + Iz - Ix)*2.5), 
                          sqrt((Iz + Ix - Iy)*2.5), 
                          sqrt((Ix + Iy - Iz)*2.5)};
    
    double gds = fmax(hd[0], fmax(hd[1],hd[2]))*2. / (double)(RES - MARGIN);
    Tuple3i res(RES, RES, RES);
    for(int i = 0;i < 3;++ i)
        while ( (res[i]/2)*gds > hd[i]*2+MARGIN*gds && res[i] >= (1<<NLEVEL) )
            res[i] /= 2;
    cout << "-------------------------------------" << endl;
    cout << " RES: " << res.x << ' ' << res.y << ' ' << res.z << endl;
    cout << " GDS: " << gds << endl;
    cout << "-------------------------------------" << endl;

    Point3d mincorer(-res.x*gds*0.5, -res.y*gds*0.5, -res.z*gds*0.5);
    EllipsoidDistFunc df(hd[0], hd[1], hd[2]);
    pstuffer = new TIsoStuffer(res, NLEVEL, gds, mincorer, 
            ALPHA_LONG, ALPHA_SHORT, df);
    pmesh = new TetMesh<double>;
    pstuffer->create_mesh(pmesh);
    TetMeshWriter_Double::write_mesh(argv[4], *pmesh);

    clean_up();
    return 0;
}

