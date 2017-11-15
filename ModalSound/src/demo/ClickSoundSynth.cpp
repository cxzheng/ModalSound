#include <limits>
#include <QApplication>
#include <QDir>
#include <QFileInfo>
#include <QKeyEvent>
#include <QAudio>
#include <QSettings>
#include <QAudioDeviceInfo>
#include <QGLViewer/qglviewer.h>
#if defined(__APPLE__) || defined(MACOSX)
#   include <glut.h>
#elif defined(__linux)
#   include <GL/glut.h>
#else
#   error ERROR Unsupported OS
#endif
#include "io/TglMeshReader.hpp"
#include "utils/macros.h"
#include "utils/term_msg.h"
#include "AudioProducer.h"

using namespace std;

class ClickSoundViewer : public QGLViewer
{
    public:
        ClickSoundViewer(const char* inifile);

        ~ClickSoundViewer()
        {
            delete audio_;
        }

    protected:
        void init();
        void draw();
        void keyPressEvent(QKeyEvent* e);
        void drawWithNames();
        void postSelection(const QPoint&);

    private:
        void init_gl();
        void load_mesh(const QString& meshfile);

        void draw_mesh();
        void draw_obj();
        void draw_triangle_normal();

    private:
        QDir                dataDir_;
        bool                wireframe_;
        int                 selTriId_;

        AudioProducer*      audio_;

        TriangleMesh<double>    mesh_;
};

// ------------------------------------------------------------

ClickSoundViewer::ClickSoundViewer(const char* inifile) : 
        wireframe_(false), selTriId_(-1), audio_(NULL)
{
    QString iniPath(inifile);
    QFileInfo checkConfig(iniPath);

    if ( !checkConfig.exists() )
    {
        PRINT_ERROR("The INI file %s doesn't exist\n", inifile);
        SHOULD_NEVER_HAPPEN(-1);
    }

    // load the settings
    QSettings settings(iniPath, QSettings::IniFormat);
    //cout << "HELLO: " << settings.value("mesh/surface_mesh").toString().toStdString() << endl;
    //cout << "       \"" << settings.value("audio/device2").toString().toStdString() << '"' << endl;

    dataDir_ = checkConfig.absoluteDir(); //.absolutePath().toStdString() << endl;

    ////// load the files
    // - triangle mesh 
    load_mesh(settings.value("mesh/surface_mesh").toString());

    audio_ = new AudioProducer(settings, dataDir_);
}

void ClickSoundViewer::load_mesh(const QString& meshfile)
{
    QString filename;
    {
        QFileInfo fInfo(meshfile);
        filename = fInfo.isRelative() ? dataDir_.filePath(meshfile) : meshfile;
    }

    PRINT_MSG("Load the mesh file: %s\n", filename.toStdString().c_str());
    if ( _FAILED( MeshObjReader::read(filename.toStdString().c_str(), mesh_) ) )
    {
        PRINT_ERROR("Cannot load the mesh file: %s\n", 
                    filename.toStdString().c_str());
        SHOULD_NEVER_HAPPEN(-1);
    }
    if ( !mesh_.has_normals() ) mesh_.generate_pseudo_normals();
}

void ClickSoundViewer::init()
{
    init_gl();      // initialize OpenGL

    resize(1024, 756);
    camera()->setZNearCoefficient(0.0001f);
    camera()->setZClippingCoefficient(100.f);
}

void ClickSoundViewer::draw()
{
    draw_mesh();
    if ( selTriId_ >= 0 ) draw_triangle_normal();
}

void ClickSoundViewer::draw_obj()
{
    const vector<Point3d>&  vtx = mesh_.vertices();
    const vector<Tuple3ui>& tgl = mesh_.surface_indices();
    const vector<Vector3d>& nml = mesh_.normals();

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)(vtx.data()));
    glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)(nml.data()));
    glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT, (const GLvoid*)(tgl.data()));
}

void ClickSoundViewer::draw_mesh()
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3f(1.f, 1.f, 0.f);
    glEnable(GL_LIGHTING);

    draw_obj();
    if ( wireframe_ )
    {
        glLineWidth(2.);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glColor3f(0.8f, 0.8f, 0.8f);
        glDisable(GL_LIGHTING);
        draw_obj();
    }
}

void ClickSoundViewer::draw_triangle_normal()
{
    const int tid = selTriId_;
    const vector<Point3d>&  vtx = mesh_.vertices();
    const vector<Tuple3ui>& tgl = mesh_.surface_indices();

    const Point3d ctr  = (vtx[tgl[tid].x] + vtx[tgl[tid].y] + vtx[tgl[tid].z])*(1./3.);
    Vector3d nml = Triangle<double>::weighted_normal( vtx[tgl[tid].x], vtx[tgl[tid].y], vtx[tgl[tid].z] );
    double area = nml.normalize2();
    const double len = 0.05;     ///// HARD coded value here!!!
    const Point3d end  = ctr + nml * len;

    glColor3f(1.f, 0.f, 0.f);
    drawArrow( qglviewer::Vec(end.x, end.y, end.z), 
               qglviewer::Vec(ctr.x, ctr.y, ctr.z),
               sqrt(area)*0.8 );
}

void ClickSoundViewer::keyPressEvent(QKeyEvent* e)
{
    const Qt::KeyboardModifiers modifiers = e->modifiers();

    if ( e->key() == Qt::Key_W && modifiers == Qt::NoButton )
    {
        wireframe_ = !wireframe_;
        update();
    }
    else
        QGLViewer::keyPressEvent(e);
}

void ClickSoundViewer::drawWithNames()
{
    glEnableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);

    glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)(mesh_.vertices().data()));
    const vector< Tuple3ui >& tgl = mesh_.triangles();
    for(size_t i = 0;i < tgl.size();++ i)
    {
        glPushName(i);
        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, (const GLvoid*)&tgl[i]); 
        glPopName();
    }
}

void ClickSoundViewer::postSelection(const QPoint&)
{
    selTriId_ = selectedName();
    //cout << "selected triangle ID:  " << selTriId_ << endl;

    //// Now synthesize sound and play
    {
        qglviewer::Vec cam = camera()->position();
        const Point3d camPos( cam.x, cam.y, cam.z );

        const vector<Point3d>&  vtx = mesh_.vertices();
        const vector<Tuple3ui>& tgl = mesh_.surface_indices();
        Vector3d nml = Triangle<double>::normal( 
                vtx[tgl[selTriId_].x], 
                vtx[tgl[selTriId_].y], 
                vtx[tgl[selTriId_].z] );
        nml.normalize();
        audio_->play( mesh_.triangle_ids(selTriId_), nml, camPos );
    }
    
    if ( selTriId_ >= 0 ) update();
}

void ClickSoundViewer::init_gl()
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

    // antialiasing
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

////////////////////////////////////////////////////////////////////////////////////////////
static void list_audio_devices()
{
    cout << "Detected audio devices: " << endl;
    foreach(const QAudioDeviceInfo &deviceInfo, 
            QAudioDeviceInfo::availableDevices(QAudio::AudioOutput)) 
        cout << "  + \"" << deviceInfo.deviceName().toStdString() << '"' << endl;
}

int main(int argc, char* argv[])
{
    if ( argc != 2 ) 
    {
        cerr << "Invalid arguments!" << endl;
        cerr << "Usage: " << argv[0] << " [list | .ini file]" << endl;
        return 1;
    }

    // list audio output devices
    if ( !strcmp(argv[1], "list") )
    {
        list_audio_devices();
        return 0;
    }

    QApplication app(argc, argv);
    ClickSoundViewer viewer(argv[1]);

    //viewer.setWindowTitle(title.c_str());
    viewer.show();

    return app.exec();
}

