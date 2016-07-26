#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/unordered_map.hpp>
#include "io/TglMeshReader.hpp"

using namespace std;

typedef TriangleMesh<double> TMesh;

static TMesh mesh;

static void check_manifold()
{
    const vector<Point3d>& vtx = mesh.vertices();
    const vector<Tuple3ui>& tgl= mesh.triangles();

    boost::unordered_map< int, boost::unordered_map<int, int> > edgecnt;
    double minlength = 1E+300;
    for(size_t i = 0;i < tgl.size();++ i)
    for(int j = 0;j < 3;++ j)
    {
        int v1 = tgl[i][j], v2 = tgl[i][(j+1)%3];
        if ( v1 > v2 ) swap(v1, v2);
        minlength = min(minlength, vtx[v1].distance(vtx[v2]));
        ++ edgecnt[v1][v2];
    }

    cout << "Mesh Check: [Manifold]" << endl;
    const boost::unordered_map<int, boost::unordered_map<int,int> >::iterator end = edgecnt.end();
    for(boost::unordered_map<int, boost::unordered_map<int,int> >::iterator it1 = edgecnt.begin();
            it1 != end; ++ it1)
    {
        const boost::unordered_map<int,int>::iterator end2 = it1->second.end();
        for(boost::unordered_map<int,int>::iterator it2 = it1->second.begin();
                it2 != end2;++ it2)
        {
            if ( it2->second != 2 )
            {
                cout << "    FAILED: " << it1->first << "--" << it2->first << " : " << it2->second << endl;
                cout << "    min len: " << minlength << endl;
            }
        }
    }
    cout << "    Succeed" << endl;
    cout << "    min len: " << minlength << endl;
}

int main(int argc, char* argv[])
{
    if ( argc != 2 ) 
    {
        cerr << "ERROR: invalid argument!" << endl;
        return 1;
    }
    MeshObjReader::read(argv[1], mesh);
    check_manifold();
    return 0;
}
