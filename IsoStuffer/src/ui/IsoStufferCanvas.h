/******************************************************************************
 *  File: IsoStufferCanvas.h
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
#ifndef UI_ISOSTUFFER_CANVAS_H
#   define UI_ISOSTUFFER_CANVAS_H

#include <QGLViewer/qglviewer.h>

class IsoStufferFrame;

class IsoStufferCanvas : public QGLViewer
{
    friend class IsoStufferFrame;

    Q_OBJECT

    public slots:
        void switch_wireframe();
        void switch_facet_disp();

    public:
        IsoStufferCanvas(QWidget* parent):QGLViewer(parent), 
                wireframe_(false), facet_(true), scale_(1), 
                sphereDisp_(0)
        { }

        void init();
        void draw();
        void keyPressEvent(QKeyEvent* e);

    private:
        void display_step1();
        void display_step0();

        void display_step23();
        void display_step46();

    private:
        IsoStufferFrame*    parent_;

        bool                wireframe_;
        bool                facet_;
        double              scale_;
        unsigned int        sphereDisp_;
};

#endif
