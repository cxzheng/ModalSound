/*
 * =====================================================================================
 *
 *       Filename:  TetViewerFrame.h
 *
 *        Version:  1.0
 *        Created:  11/17/10 16:38:15
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef  TET_VIEWER_FRAME_INC
#define  TET_VIEWER_FRAME_INC

#include <QApplication>
#include <QLabel>
#if QT_VERSION >= 0x040000
#include "ui_tetviewer.h"
#else
#error Only Qt with the version higher than 4.0 is supported right now
#endif
#include "ModesVisParamsDialog.h"
#include "geometry/FixedVtxTetMesh.hpp"

class TetViewerFrame : public QMainWindow, private Ui_TetViewerFrame
{
    friend class TetViewerCanvas;

    Q_OBJECT

    public slots:
        void open();
        void export_bin_tet();
        void export_abaqus_tet();
        void check_useless_vtx();
        void load_modes();
        void modes_params();

        void set_mode_sel(int mId);
        void mode_scale_changed(double s);

        void mode_vis_finished(int)
        {   canvas->set_disp_mode(TetViewerCanvas::SHOW_MESH); }

    public:
        typedef FixedVtxTetMesh<double>     TMesh;

        TetViewerFrame():modesVisDiag_(this), mesh_(NULL), modeSel_(0)
        {
            setupUi(this);
            canvas->parent_ = this;

            QObject::connect(actionOpen, SIGNAL(triggered()), this, SLOT(open()));
            QObject::connect(actionWireframe, SIGNAL(toggled(bool)), canvas, SLOT(toggle_wireframe(bool)));
            QObject::connect(actionMeshInfo, SIGNAL(toggled(bool)), canvas, SLOT(toggle_show_info(bool)));
            QObject::connect(actionBinaryTetFormat, SIGNAL(triggered()), this, SLOT(export_bin_tet()));
            QObject::connect(actionAbaqusTetFormat, SIGNAL(triggered()), this, SLOT(export_abaqus_tet()));
            QObject::connect(actionCheckUselessVertex, SIGNAL(triggered()), this, SLOT(check_useless_vtx()));

            QObject::connect(actionLoadModes, SIGNAL(triggered()), this, SLOT(load_modes()));
            QObject::connect(actionModeShapes, SIGNAL(triggered()), this, SLOT(modes_params()));
            QObject::connect(&modesVisDiag_, SIGNAL(scaleChanged(double)), this, SLOT(mode_scale_changed(double)));
            QObject::connect(&modesVisDiag_, SIGNAL(finished(int)), this, SLOT(mode_vis_finished(int)));

            statusbar->addWidget(new QLabel("  Press \"H\" for help   ", statusbar));
        }
        ~TetViewerFrame()
        {   delete mesh_; }

    private:
        void update_mesh(TMesh* msh);
        void update_normals();

    private:
        ModesVisParamsDialog    modesVisDiag_;

        TMesh*                  mesh_;
        std::vector<Point3d>    vtx_;       // vertices in the tet mesh
        std::vector<Vector3d>   nml_;       // vertex normal. if the vertex is not on the surface
                                            // its normal is undefined.

        /* --------- for mode visualization -------- */
        int                                 modeSel_;
        double                              modeScale_;
        std::vector< std::vector<double> >  modes_;
        std::vector< double >               maxMAmp_;
};

#endif   /* ----- #ifndef TET_VIEWER_FRAME_INC  ----- */

