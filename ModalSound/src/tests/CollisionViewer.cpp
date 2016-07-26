#include <iostream>
#include <QApplication>
#include <QKeyEvent>
#include <QGLViewer/qglviewer.h>
#if defined(__APPLE__) || defined(MACOSX)
#   include <glut.h>
#elif defined(__linux)
#   include <GL/glut.h>
#else
#   error ERROR Unsupported OS
#endif

#include "generic/null_type.hpp"
#include "io/TglMeshReader.hpp"
#include "collision/RigidCollisionDetector.h"

using namespace std;

class CollisionViewer : public QGLViewer
{
    public:
        CollisionViewer(const TriangleMesh<double>* meshA,
                        const TriangleMesh<double>* meshB)
        { 
            mesh_[0] = meshA;
            mesh_[1] = meshB;
            RigidCollisionDetector detector(meshA, meshB);
            int ret = detector.detect( intersectTriPairs_ );
            printf("%d intersecting tri. pairs are detected.\n", ret);
        }

    protected:
        void init();
        void draw();
        //void keyPressEvent(QKeyEvent* e);

    private:
        void init_gl();
        void draw_mesh(int id) const;

    private:
        const TriangleMesh<double>* mesh_[2];
        vector< pair<int, int> >    intersectTriPairs_;

};

void CollisionViewer::draw_mesh(int id) const
{
    const vector<Tuple3ui>& tgl = mesh_[id]->surface_indices();
    const vector<Point3d>&  vtx = mesh_[id]->vertices();
    const vector<Vector3d>& nml = mesh_[id]->normals();

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)(&(vtx[0])));
    glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)(&(nml[0])));

    glDisable(GL_LIGHTING);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    //glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT, (const GLvoid*)&tgl[0]);

    glEnable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    for(size_t i = 0;i < intersectTriPairs_.size();++ i)
    {
        const int tid = id == 0 ?  intersectTriPairs_[i].first : intersectTriPairs_[i].second;
        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, (const GLvoid*)&tgl[tid]);
    }

#if 1
    glDisable(GL_LIGHTING);
    glColor3f(1.f, 1.0f, 1.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    for(size_t i = 0;i < intersectTriPairs_.size();++ i)
    {
        const int tid = id == 0 ?  intersectTriPairs_[i].first : intersectTriPairs_[i].second;
        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, (const GLvoid*)&tgl[tid]);
    }
#endif
}

void CollisionViewer::init_gl()
{
#ifdef __linux
    int dummy = 0;
    glutInit(&dummy, NULL);
#endif

    const GLfloat GLOBAL_AMBIENT[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    const GLfloat SPECULAR_COLOR[] = { 0.6f, 0.6f, 0.6f, 1.0 };

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, GLOBAL_AMBIENT);
    const GLfloat ambientLight[] = { 0.f, 0.f, 0.f, 1.0f };
    const GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8f, 1.0f };
    const GLfloat specularLight[] = { 1.f, 1.f, 1.f, 1.0f };
    const GLfloat position[] = { -0.5f, 1.0f, 0.4f, 1.0f };

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, SPECULAR_COLOR);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    // anti-aliasing
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void CollisionViewer::init()
{
    init_gl();

    resize(1024, 756);
    camera()->setZNearCoefficient(0.0001f);
    camera()->setZClippingCoefficient(100.f);
}

void CollisionViewer::draw()
{
    glLineWidth(2.);
    glColor3f(0.8f, 0.0f, 0.0f);
    draw_mesh(0);

    glColor3f(0.0f, 0.0f, 0.8f);
    draw_mesh(1);
}
// ------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    if ( argc != 3 ) 
    {
        cerr << "Invalid arguments!" << endl;
        cerr << "Usage: " << argv[0] << " [first mesh .obj file] [second mesh .obj file]" << endl;
        return 1;
    }

    TriangleMesh<double> meshA, meshB;

    //// load mesh
    MeshObjReader::read(argv[1], meshA);
    if ( !meshA.has_normals() ) meshA.generate_pseudo_normals();
    MeshObjReader::read(argv[2], meshB);
    if ( !meshB.has_normals() ) meshB.generate_pseudo_normals();

    QApplication app(argc, argv);

    CollisionViewer viewer(&meshA, &meshB);

    string title = "Mesh File: " + string(argv[1]) + string("/") + string(argv[2]);
    viewer.setWindowTitle(title.c_str());
    viewer.show();

    return app.exec();
}
