/******************************************************************************
 *  File: IsoStufferFrame.h
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
#ifndef UI_ISOSTUFFER_FRAME
#   define UI_ISOSTUFFER_FRAME

#include <QApplication>
#include <QLabel>

#if QT_VERSION >= 0x040000
#include "ui_IsoStufferFrame.h"
#else
#error Only Qt with the version higher than 4.0 is supported right now
#endif
#include "IsoStuffer.hpp"
#include "geometry/Point3.hpp"
#include "ParamsDialog.h"
#include "DistValProc.h"

struct SphereDistFunc
{
    double operator() (const Point3d& pt)
    {
        //return pt.norm() - 1.00123456789;
        return pt.x*pt.x*2 + pt.y*pt.y*1.5 + pt.z*pt.z*1.0 - 1.0;
    }
};

class IsoStufferCanvas;

class IsoStufferFrame : public QMainWindow, private Ui_IsoStufferFrame
{
    friend class IsoStufferCanvas;

    Q_OBJECT

    public slots:
        void next_step();
        void show_params() {   params_.show(); }
        void load_mesh();
        void save_mesh();
        void export_to_stellar();
        void convert_to_tetmesh();
        void check_mesh();

        void reload_mesh()
        {
#if QT_VERSION >= 0x050000
            load_tgl_mesh(curFile_.toLatin1().data());
#else
            load_tgl_mesh(curFile_.toAscii().data());
#endif
        }

    protected:
        void closeEvent(QCloseEvent *event);

    public:
        IsoStufferFrame():curStep_(0), pstuffer_(NULL), 
                ptglmesh_(NULL), pmesh_(NULL), dvproc_(NULL),
                pdistfunc_(NULL)
        {
            setupUi(this);
            canvas->parent_ = this;

            QObject::connect(actionStep, SIGNAL(triggered()), this, SLOT(next_step()));
            QObject::connect(actionWireframe, SIGNAL(triggered()), canvas, SLOT(switch_wireframe()));
            QObject::connect(actionFacet, SIGNAL(triggered()), canvas, SLOT(switch_facet_disp()));
            QObject::connect(actionParams, SIGNAL(triggered()), this, SLOT(show_params()));
            QObject::connect(actionLoad, SIGNAL(triggered()), this, SLOT(load_mesh()));
            QObject::connect(actionReload, SIGNAL(triggered()), this, SLOT(reload_mesh()));
            QObject::connect(actionConvert, SIGNAL(triggered()), this, SLOT(convert_to_tetmesh()));
            QObject::connect(actionSaveMesh, SIGNAL(triggered()), this, SLOT(save_mesh()));
            QObject::connect(actionExportStellar, SIGNAL(triggered()), this, SLOT(export_to_stellar()));
            QObject::connect(actionCheckMesh, SIGNAL(triggered()), this, SLOT(check_mesh()));
            statusbar->addWidget(new QLabel("  Press \"H\" for help   ", statusbar));            

            //// for testing
            //pstuffer_ = new TIsoStuffer(Tuple3i(32, 32, 32), 4, 
            //        3./32., Point3d(-1.5, -1.5, -1.5), 0.25, 0.41189,
            //        SphereDistFunc());
            //pstuffer_ = new TIsoStuffer(Tuple3i(64, 64, 64), 5, 
            //        2.4/64., Point3d(-1.2, -1.2, -1.2), 0.25, 0.41189,
            //        SphereDistFunc());
        }
        
        ~IsoStufferFrame()
        {   reset(); }

    private:
        void reset();
        void load_tgl_mesh(const char* file);

    private:
        typedef DistFunc<double, DistValProc>           TDistFunc;
        typedef IsoStuffing::IsoStuffer<TDistFunc>      TIsoStuffer;

        int                     curStep_;
        TIsoStuffer*            pstuffer_;
        TriangleMesh<double>*   ptglmesh_;
        TetMesh<double>*        pmesh_;
        DistValProc*            dvproc_;
        TDistFunc*              pdistfunc_;
        ParamsDialog            params_;

        QString                 curFile_;
};

#endif
