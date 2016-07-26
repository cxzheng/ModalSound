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

#include "geometry/OBBTree.hpp"
#include "generic/null_type.hpp"
#include "io/TglMeshReader.hpp"

using namespace std;

typedef OBBTree< double, TriangleMesh<double> > TOBBTree;

static TriangleMesh<double> mesh;
static TOBBTree*            tree = NULL;

class OBBViewer : public QGLViewer
{
    public:
        OBBViewer(const TriangleMesh<double>* mesh, TOBBTree* tree):
                mesh_(mesh), tree_(tree), level_(0), 
                wireframe_(false), showMesh_(true) { }

    protected:
        void init();
        void draw();
        void keyPressEvent(QKeyEvent* e);

    private:
        void init_gl();
        void draw_mesh() const;
        void draw_bounding_box(const TOBBTree::TNode* node) const;
        void toggle_wireframe(bool wf);

    private:
        const TriangleMesh<double>* mesh_;
        const TOBBTree*             tree_;

        int                         level_;
        bool                        wireframe_;
        bool                        showMesh_;
};

void OBBViewer::draw_mesh() const
{
    const vector<Tuple3ui>& tgl = mesh_->surface_indices();
    const vector<Point3d>&  vtx = mesh_->vertices();
    const vector<Vector3d>& nml = mesh_->normals();

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)(&(vtx[0])));
    glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)(&(nml[0])));
    glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT, (const GLvoid*)&tgl[0]);
}

void OBBViewer::init_gl()
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

void OBBViewer::init()
{
    init_gl();

    resize(1024, 756);
    camera()->setZNearCoefficient(0.0001f);
    camera()->setZClippingCoefficient(100.f);
}

void OBBViewer::draw()
{
    if ( showMesh_ )
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f(0.f, 0.5f, 1.f);
        glEnable(GL_LIGHTING);
        draw_mesh();

        if ( wireframe_ )
        {
            glLineWidth(2.);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glColor3f(0.8f, 0.8f, 0.8f);
            glDisable(GL_LIGHTING);
            draw_mesh();
        }
    }
    draw_bounding_box(tree_->root());
}

static const int _DIR[][2] = { {-1, -1}, {-1, 1}, {1, 1}, {1, -1}};
void OBBViewer::draw_bounding_box(const TOBBTree::TNode* node) const
{
    if ( node->level() == level_ )
    {
        //// draw current tree node
        const Matrix3<double>& R = node->R();                                                               
        const Tuple3<double>&  r = node->r();                                                               
        const Point3<double>&  c = node->c();                                                               

        glColor3f(1.f, 0.f, 0.f);
        glBegin(GL_LINES);

        for(int i = 0;i < 4;++ i)
        {
            Vector3<double> dd1(_DIR[i][0]*r[0], _DIR[i][1]*r[1], -r[2]);
            Vector3<double> dd2(_DIR[i][0]*r[0], _DIR[i][1]*r[1],  r[2]);                                   
            glVertex3dv(c + R*dd1);
            glVertex3dv(c + R*dd2); 
    
            dd2.set(_DIR[(i+1)%4][0]*r[0], _DIR[(i+1)%4][1]*r[1], -r[2]);                                 

            glVertex3dv(c + R*dd1);
            glVertex3dv(c + R*dd2);

            dd1[2] *= -1;
            dd2[2] *= -1;                                                                                 
            glVertex3dv(c + R*dd1);                                                                       
            glVertex3dv(c + R*dd2);
        }
    
        glEnd();
        return;
    }

    const TOBBTree::TNode** c = node->children();
    if ( c[0] ) draw_bounding_box(c[0]);
    if ( c[1] ) draw_bounding_box(c[1]);
}

void OBBViewer::keyPressEvent(QKeyEvent* e)
{
    const Qt::KeyboardModifiers modifiers = e->modifiers();

    if ( e->key() == Qt::Key_W && modifiers == Qt::NoButton )
        toggle_wireframe(!wireframe_);
    else if ( e->key() == Qt::Key_Less )
    {   
        -- level_;                                                                                       
        level_ = std::max(0, level_);
        updateGL();
    }   
    else if ( e->key() == Qt::Key_Greater )
    {   
        ++ level_;
        updateGL();
    }   
    else
        QGLViewer::keyPressEvent(e);
}

void OBBViewer::toggle_wireframe(bool wf)
{
    wireframe_ = wf;
    updateGL();
}

// ------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    if ( argc != 2 ) 
    {
        cerr << "Invalid arguments!" << endl;
        cerr << "Usage: " << argv[0] << " [mesh .obj file]" << endl;
        return 1;
    }

    //// load mesh
    MeshObjReader::read(argv[1], mesh);
    if ( !mesh.has_normals() ) mesh.generate_pseudo_normals();

    tree = new TOBBTree(&mesh);

    QApplication app(argc, argv);
    OBBViewer viewer(&mesh, tree);

    string title = "Mesh File: " + string(argv[1]);
    viewer.setWindowTitle(title.c_str());
    viewer.show();

    return app.exec();
}
