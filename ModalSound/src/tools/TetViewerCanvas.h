/*
 * =====================================================================================
 *
 *       Filename:  TetViewerCanvas.h
 *
 *        Version:  1.0
 *        Created:  11/17/10 16:41:25
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef  TET_VIEWER_CANVAS_INC
#define  TET_VIEWER_CANVAS_INC

#include <QGLViewer/qglviewer.h>
#include "colormap/ColorMap.h"

class TetViewerFrame;

class TetViewerCanvas : public QGLViewer
{
    friend class TetViewerFrame;

    Q_OBJECT

    public slots:
        void toggle_wireframe(bool wf);
        void toggle_show_info(bool si);

    public:
        TetViewerCanvas(QWidget* parent):QGLViewer(parent), 
                wireframe_(false), showMshInfo_(false)
        { }

    protected:
        void init();
        void draw();
        void keyPressEvent(QKeyEvent* e);

    private:
        enum DISP_MODE
        {
            SHOW_MESH = 0,
            SHOW_MODES
        };

        void init_gl();
        void draw_obj();
        void show_mesh_info();

        void draw_mesh();
        void draw_modes();
        void set_disp_mode(DISP_MODE d)
        {   dispMode_ = d; }

        void dump_mode_colormap(const char* fname);

    private:
        bool            wireframe_;
        bool            showMshInfo_;
        DISP_MODE       dispMode_;
        Um4ColorMap     cMap_;

        TetViewerFrame* parent_;
};

#endif   /* ----- #ifndef TET_VIEWER_CANVAS_INC  ----- */
