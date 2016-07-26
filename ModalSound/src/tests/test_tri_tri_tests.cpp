#include <iostream>
#include <math.h>
#include "geometry/Point3.hpp"
#include "geometry/tritri.h"

using namespace std;

static void test1()
{
    Point3d p0(0., 0., 0.);
    Point3d p1(1., 0., 0.);
    Point3d p2(0., 0., 1.);

    Point3d q0(0.1,  1., 0.1);
    Point3d q1(0.1, -1., 0.1);
    Point3d q2(sqrt(2.)+0.1, -1., sqrt(2.)+0.1);

    Point3d source, target;
    int coplanar;

    int ret = tri_tri_intersection_test_3d(
            (const double*)&p0, (const double*)&p1, (const double*)&p2,
            (const double*)&q0, (const double*)&q1, (const double*)&q2,
            &coplanar,
            (double*)&source, (double*)&target);
    cout << "RET: " << ret << endl;
    cout << "CP:  " << coplanar << endl;
    cout << "SRC: " << source << endl;
    cout << "TGT: " << target << endl;
    cout << "     " << tri_tri_overlap_test_3d(
            (const double*)&p0, (const double*)&p1, (const double*)&p2,
            (const double*)&q0, (const double*)&q1, (const double*)&q2)
         << endl;
}

static void test2()
{
    Point3d p0(0., 0., 0.);
    Point3d p1(1., 0., 0.);
    Point3d p2(0., 0., 1.);

    Point3d q0(0.2, 0., 0.2);
    Point3d q1(1.2, 0., 0.2);
    Point3d q2(0.2, 0., 1.2);

    Point3d source, target;
    int coplanar;

    int ret = tri_tri_intersection_test_3d(
            (const double*)&p0, (const double*)&p1, (const double*)&p2,
            (const double*)&q0, (const double*)&q1, (const double*)&q2,
            &coplanar,
            (double*)&source, (double*)&target);
    cout << "RET: " << ret << endl;
    cout << "CP:  " << coplanar << endl;
    cout << "SRC: " << source << endl;
    cout << "TGT: " << target << endl;
    cout << "     " << tri_tri_overlap_test_3d(
            (const double*)&p0, (const double*)&p1, (const double*)&p2,
            (const double*)&q0, (const double*)&q1, (const double*)&q2)
         << endl;
}
int main(int argc, char* argv[])
{
    test1();
    test2();
    return 0;
}
