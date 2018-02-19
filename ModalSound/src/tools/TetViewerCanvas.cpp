/*
 * =====================================================================================
 *
 *       Filename:  TetViewerCanvas.cpp
 *
 *        Version:  1.0
 *        Created:  11/17/10 17:13:17
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#include "TetViewerCanvas.h"

#if defined(__APPLE__) || defined(MACOSX)
#   include <glut.h>
#elif defined(__linux)
#   include <GL/glut.h>
#else
#   error ERROR Unsupported OS
#endif

#include <QKeyEvent>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TetViewerFrame.h"
#include "geometry/TetMesh.hpp"

using namespace std;

void TetViewerCanvas::init()
{
    init_gl();

    setKeyDescription(Qt::Key_W, "Toggles wire frame display");
    setKeyDescription(Qt::Key_I, "Toggles simulation information");
    setKeyDescription(Qt::Key_I, "Advance single step");

    camera()->setZNearCoefficient(0.0001f);
    camera()->setZClippingCoefficient(100.f);

    //setBackgroundColor(QColor::fromRgbF(1.f,1.f,1.f));
}

void TetViewerCanvas::init_gl()
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

void TetViewerCanvas::toggle_wireframe(bool wf)
{
    wireframe_ = wf;
    update();
}
void TetViewerCanvas::toggle_show_info(bool si)
{
    showMshInfo_ = si;
    update();
}

void TetViewerCanvas::dump_mode_colormap(const char* fname)
{
    ofstream fout(fname);
    if ( fout.fail() ) 
    {
        cerr << "ERROR: Fail to open the file: " << fname << endl;
        return;
    }

    const vector< Point3d >&  vtx = parent_->mesh_->vertices();
    const vector< Tuple3ui >& tgl = parent_->mesh_->surface_indices();
    const int nfv = parent_->mesh_->num_fixed_vertices();

    cMap_.set_interpolation_range(0., parent_->maxMAmp_[parent_->modeSel_]*parent_->modeScale_);

    fout << "ply" << endl;
    fout << "format ascii 1.0" << endl;
    fout << "element vertex " << vtx.size() << endl;
    fout << "property float x" << endl;
    fout << "property float y" << endl;
    fout << "property float z" << endl;
    fout << "property uchar red" << endl;
    fout << "property uchar green" << endl;
    fout << "property uchar blue" << endl;
    fout << "element face " << tgl.size() << endl;
    fout << "property list uchar int vertex_indices" << endl;
    fout << "end_header" << endl;

    fout << setprecision(8);
    for(size_t i = 0;i < vtx.size();++ i)
    {
        fout << vtx[i].x << ' ' << vtx[i].y << ' ' << vtx[i].z; // << endl;
        if ( (int)i >= nfv )
        {
            const double* ptr = &(parent_->modes_[parent_->modeSel_][(i-nfv)*3]);
            Tuple3f c = cMap_.get_interpolated_color(parent_->modeScale_*
                    sqrt(M_SQR(ptr[0])+M_SQR(ptr[1])+M_SQR(ptr[2])));
            fout << ' ' << fmin(int(c.x*255.), 255)
                 << ' ' << fmin(int(c.y*255.), 255)
                 << ' ' << fmin(int(c.z*255.), 255) 
                 << endl;
        }
        else
            fout << "3 255 255 255" << endl;
    }

    for(size_t i = 0;i < tgl.size();++ i)
    {
        fout << "3 " << tgl[i].x << ' ' << tgl[i].y << ' ' << tgl[i].z << endl;
    }
    fout.close();
}

void TetViewerCanvas::keyPressEvent(QKeyEvent* e)
{
    const Qt::KeyboardModifiers modifiers = e->modifiers();
    
    if ( e->key() == Qt::Key_W && modifiers == Qt::NoButton )
        toggle_wireframe(!wireframe_);
    else if ( e->key() == Qt::Key_I && modifiers == Qt::NoButton )
    {
        printf("=========== Mesh Statistics ===========\n");
        printf(" # of vertices:          %d\n", (int)parent_->mesh_->num_vertices());
        printf(" # of tetrahedron:       %d\n", (int)parent_->mesh_->num_tets());
        printf(" # of surface triangles: %d\n", (int)parent_->mesh_->num_surface_tgls());
        printf(" # of free vertices:     %d\n", (int)parent_->mesh_->num_free_vertices());
        printf(" # of fixed vertices:    %d\n", (int)parent_->mesh_->num_fixed_vertices());
        printf("=======================================\n");
    }
    else if ( e->key() == Qt::Key_D && modifiers == Qt::NoButton )
    {
        if ( dispMode_ == SHOW_MODES ) 
        {
            dump_mode_colormap("cmap.ply");
        }
    }
    else
        QGLViewer::keyPressEvent(e);
}

void TetViewerCanvas::show_mesh_info()
{
    glColor3f(1.f, 1.f, 1.f);
    drawText(width()-170, 20, QString("# of vtx:       %1").arg(
                parent_->mesh_->num_vertices(), 0));
    drawText(width()-170, 40, QString("# of tets:      %1").arg(
                parent_->mesh_->num_tets(), 0));
    drawText(width()-170, 60, QString("# of surf tgl:  %1").arg(
                parent_->mesh_->num_surface_tgls(), 0));
    drawText(width()-170, 80, QString("# of free vtx:  %1").arg(
                parent_->mesh_->num_free_vertices(), 0));
    drawText(width()-170, 100, QString("# of fixed vtx: %1").arg(
                parent_->mesh_->num_fixed_vertices(), 0));
    if ( !parent_->maxMAmp_.empty() )
    drawText(width()-170, 120, QString("# of vib. modes: %1").arg(
                parent_->maxMAmp_.size(), 0));
}

void TetViewerCanvas::draw_mesh()
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3f(0.f, 1.f, 0.f);
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

void TetViewerCanvas::draw_modes()
{
    const vector< Point3d >&  vtx = parent_->mesh_->vertices();
    const vector< Tuple3ui >& tgl = parent_->mesh_->surface_indices();
    const int nfv = parent_->mesh_->num_fixed_vertices();

    cMap_.set_interpolation_range(0., parent_->maxMAmp_[parent_->modeSel_]*parent_->modeScale_);

    vector< Tuple3f > colors(vtx.size());
    vector< Point3d > vtxNow(vtx.size());

    for(size_t i = 0;i < vtx.size();++ i)
    {
        vtxNow[i] = vtx[i];
        if ( (int)i >= nfv )
        {
            const double* ptr = &(parent_->modes_[parent_->modeSel_][(i-nfv)*3]);

            vtxNow[i].x += parent_->modeScale_*ptr[0];
            vtxNow[i].y += parent_->modeScale_*ptr[1];
            vtxNow[i].z += parent_->modeScale_*ptr[2];

            colors[i] = cMap_.get_interpolated_color(parent_->modeScale_*
                    sqrt(M_SQR(ptr[0])+M_SQR(ptr[1])+M_SQR(ptr[2])));
        }
        else
            colors[i].set(1.f, 1.f, 1.f);
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)(vtxNow.data()));
    glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)(parent_->nml_.data()));
    glColorPointer(3, GL_FLOAT, 0, (const GLvoid*)(colors.data()));

    // draw it
    glEnable(GL_LIGHTING);
    glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT, (const GLvoid*)tgl.data());

    glDisableClientState(GL_COLOR_ARRAY);
}

void TetViewerCanvas::draw()
{
    if ( !parent_->mesh_ ) return;

    switch ( dispMode_ )
    {
        case SHOW_MESH:
            draw_mesh();
            break;
        case SHOW_MODES:
            draw_modes();
            break;
    }

    if ( showMshInfo_ ) show_mesh_info();
}

void TetViewerCanvas::draw_obj()
{
    const vector<Tuple3ui>& tgl = parent_->mesh_->surface_indices();
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)(&(parent_->vtx_[0])));
    glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)(&(parent_->nml_[0])));
    glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT, (const GLvoid*)&tgl[0]);
}

