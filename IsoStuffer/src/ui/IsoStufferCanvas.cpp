/******************************************************************************
 *  File: IsoStufferCanvas.cpp
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
#include "IsoStufferCanvas.h"
#if defined(__APPLE__) || defined(MACOSX)
#   include <glut.h>
#elif defined(__linux)
#   include <GL/glut.h>
#else
#   error ERROR Unsupported OS
#endif

#include "IsoStufferFrame.h"
#include <QKeyEvent>

using namespace std;

void IsoStufferCanvas::init()
{
#ifdef __linux
    int dummy = 0;
    glutInit(&dummy, NULL);
#endif

    camera()->setZNearCoefficient(0.0001f);
    camera()->setZClippingCoefficient(100.f);

    setKeyDescription(Qt::Key_W, "Turn on/off wire frame display");
    glShadeModel(GL_SMOOTH);
}

void IsoStufferCanvas::keyPressEvent(QKeyEvent* e)
{
    const Qt::KeyboardModifiers modifiers = e->modifiers();
    
    if ( e->key() == Qt::Key_W && modifiers == Qt::NoButton )
        switch_wireframe();
    if ( e->key() == Qt::Key_F && modifiers == Qt::NoButton )
        switch_facet_disp();
    else
        QGLViewer::keyPressEvent(e);
}

void IsoStufferCanvas::switch_facet_disp()
{
    facet_ = !facet_;
    update();
}

void IsoStufferCanvas::switch_wireframe()
{
    wireframe_ = !wireframe_;
    if ( wireframe_ )
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT, GL_FILL);
    update();
}

void IsoStufferCanvas::draw()
{
    switch ( parent_->curStep_ )
    {
        case 0:
            display_step0();
            break;

        case 1:
            display_step1();
            break;
        case 2:
        case 3:
            display_step23();
            break;
        case 4:
        case 6:
            display_step46();
            break;
    };
}

void IsoStufferCanvas::display_step1()
{
    if ( !sphereDisp_ )
    {
        sphereDisp_ = glGenLists(1);
        glNewList(sphereDisp_, GL_COMPILE);
        glutSolidSphere(1, 20, 20);
        glEndList();
    }

    glEnable(GL_LIGHTING);
    glColor3f(0, 1, 1);

    const IsoStuffing::OctTree& tree = parent_->pstuffer_->oct_tree();
    const Tuple3i& res  = tree.resolution(0);
    const IsoStuffing::OctTree::TOctLevel& level = tree.level(0);
    for(int iz = 0;iz < res.z;++ iz)
    for(int iy = 0;iy < res.y;++ iy)
    for(int ix = 0;ix < res.x;++ ix)
        if ( level[iz][iy][ix] )
        {
            Point3d pts(ix, iy, iz);
            pts += 0.5;
            pts *= parent_->pstuffer_->smallest_grid_size();
            pts += parent_->pstuffer_->min_corner();
            glPushMatrix();
            glTranslated(pts.x, pts.y, pts.z);
            glScaled(scale_, scale_, scale_);
            glCallList(sphereDisp_);
            glPopMatrix();
        }
}

void IsoStufferCanvas::display_step23()
{
    if ( !sphereDisp_ )
    {
        sphereDisp_ = glGenLists(1);
        glNewList(sphereDisp_, GL_COMPILE);
        glutSolidSphere(0.003, 20, 20);
        glEndList();
    }

    glEnable(GL_LIGHTING);

    const IsoStuffing::OctTree& tree = parent_->pstuffer_->oct_tree();
    const vector<Tuple3i>& ress = tree.resolutions();
    const double cdelta = 1. / (double)ress.size();

    double sz = parent_->pstuffer_->smallest_grid_size();
    for(size_t i = 0;i < ress.size();++ i, sz *= 2.)
    {
        const IsoStuffing::OctTree::TOctLevel& level = tree.level((int)i);
        glColor3f(1, (float)((i+1)*cdelta), (float)((i+1)*cdelta));

        for(int iz = 0;iz < ress[i].z;++ iz)
        for(int iy = 0;iy < ress[i].y;++ iy)
        for(int ix = 0;ix < ress[i].x;++ ix)
            if ( level[iz][iy][ix] )
            {
                Point3d pts(ix, iy, iz);
                pts += 0.5;
                pts *= sz;
                pts += parent_->pstuffer_->min_corner();
                if ( pts.z < 0.4 && pts.z > -0.1 )
                //if ( i > ress.size() - 2 )
                {
                    glPushMatrix();
                    glTranslated(pts.x, pts.y, pts.z);
                    //glScalef((float)(i+1), (float)(i+1), (float)(i+1));
                    float ss = (float)(1 << i);
                    glScaled(scale_*ss, scale_*ss, scale_*ss);
                    glCallList(sphereDisp_);
                    glPopMatrix();
                }
            }
    }
}

void IsoStufferCanvas::display_step46()
{
    if ( !parent_->pmesh_ ) return;

    const vector<Point3d>&      vtx = parent_->pmesh_->vertices();
    const vector<Tuple3ui>&     tgl = parent_->pmesh_->surface_indices();

    if ( !facet_ )
    {
        const vector<Vector3d>&     nml = parent_->pmesh_->normals();

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);

        glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)&vtx[0]);
        glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)&nml[0]);

        glColor3f(0, 1, 1);
        glEnable(GL_LIGHTING);
        //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT, (const GLvoid*)&tgl[0]);
    }
    else
    {
        glColor3f(0, 1, 1);
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        for(size_t i = 0;i < tgl.size();++ i)
        {
            Vector3d nl = Triangle<double>::normal(vtx[tgl[i][0]], vtx[tgl[i][1]], vtx[tgl[i][2]]);
            nl.normalize();

            glNormal3dv(nl);
            glVertex3dv(vtx[tgl[i][0]]);
            glNormal3dv(nl);
            glVertex3dv(vtx[tgl[i][1]]);
            glNormal3dv(nl);
            glVertex3dv(vtx[tgl[i][2]]);
        }
        glEnd();
    }
}

void IsoStufferCanvas::display_step0()
{
    if ( !parent_->ptglmesh_ ) return;

    const vector<Point3d>&      vtx = parent_->ptglmesh_->vertices();
    const vector<Tuple3ui>&     tgl = parent_->ptglmesh_->surface_indices();

    if ( !facet_ )
    {
        const vector<Vector3d>&     nml = parent_->ptglmesh_->normals();

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);

        glVertexPointer(3, GL_DOUBLE, 0, (const GLvoid*)&vtx[0]);
        glNormalPointer(GL_DOUBLE, 0, (const GLvoid*)&nml[0]);

        glColor3f(0, 1, 1);
        glEnable(GL_LIGHTING);
        //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawElements(GL_TRIANGLES, tgl.size()*3, GL_UNSIGNED_INT, (const GLvoid*)&tgl[0]);
    }
    else
    {
        glColor3f(0.8, 0.8, 0.8);
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        for(size_t i = 0;i < tgl.size();++ i)
        {
            Vector3d nl = Triangle<double>::normal(vtx[tgl[i][0]], vtx[tgl[i][1]], vtx[tgl[i][2]]);
            nl.normalize();

            glNormal3dv(nl);
            glVertex3dv(vtx[tgl[i][0]]);

            glNormal3dv(nl);
            glVertex3dv(vtx[tgl[i][1]]);

            glNormal3dv(nl);
            glVertex3dv(vtx[tgl[i][2]]);
        }
        glEnd();
    }

}
