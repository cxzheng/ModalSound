/*
 * =====================================================================================
 *
 *       Filename:  TetViewerFrame.cpp
 *
 *        Version:  1.0
 *        Created:  11/17/10 16:43:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#include "TetViewerFrame.h"
#include <QFileDialog>
#include "io/StellarIO.hpp"
#include "io/TetMeshReader.hpp"
#include "io/TetMeshWriter.hpp"
#include "utils/term_msg.h"

using namespace std;

void TetViewerFrame::open()
{
    QString file = QFileDialog::getOpenFileName(this,
            "Select the mesh file", ".", 
            "mesh bin file (*.tet);;mesh txt file (*.node)");
    if ( file.isEmpty() ) return;
    
    TMesh* msh = new TMesh;
    if ( file.endsWith(".tet", Qt::CaseInsensitive) )
    {
#if QT_VERSION >= 0x050000
        FV_TetMeshLoader_Double::load_mesh(file.toLatin1().data(), *msh);
        //TetMeshLoader_Double::load_mesh(file.toLatin1().data(), *msh);
#else
        FV_TetMeshLoader_Double::load_mesh(file.toAscii().data(), *msh);
        //TetMeshLoader_Double::load_mesh(file.toAscii().data(), *msh);
#endif
    }
    else if ( file.endsWith(".node", Qt::CaseInsensitive) )
    {
#if QT_VERSION >= 0x050000
        StellarTetMeshLoader::load_mesh(file.left(file.length()-5).toLatin1().data(), msh);
#else
        StellarTetMeshLoader::load_mesh(file.left(file.length()-5).toAscii().data(), msh);
#endif
    }
    msh->init();
    msh->update_surface();

    update_mesh(msh);
}

void TetViewerFrame::export_bin_tet()
{
    if ( !mesh_ ) return;

    QString file = QFileDialog::getSaveFileName(this,
            "Binary Mesh file name", ".", "bin tet file (*.tet);;All (*)");
    if ( file.isEmpty() ) return;
#if QT_VERSION >= 0x050000
    FV_TetMeshWriter_Double::write_mesh(file.toLatin1().data(), *mesh_);
    //TetMeshWriter_Double::write_mesh(file.toLatin1().data(), *mesh_);
#else
    FV_TetMeshWriter_Double::write_mesh(file.toAscii().data(), *mesh_);
    //TetMeshWriter_Double::write_mesh(file.toAscii().data(), *mesh_);
#endif
}

void TetViewerFrame::export_abaqus_tet()
{
    if ( !mesh_ ) return;
    QString file = QFileDialog::getSaveFileName(this,
            "Abaqus Mesh file name", ".", "abaqus tet file (*.aba);;All (*)");
    if ( file.isEmpty() ) return;
#if QT_VERSION >= 0x050000
    FV_AbaqusMeshWriter::write_mesh(file.toLatin1().data(), *mesh_);
#else
    FV_AbaqusMeshWriter::write_mesh(file.toAscii().data(), *mesh_);
#endif
}

void TetViewerFrame::update_mesh(TMesh* msh)
{
    if ( mesh_ ) delete mesh_;

    vtx_.resize(msh->num_vertices());
    const vector<Point3d>& vs = msh->vertices();
    memcpy(&vtx_[0], &vs[0], sizeof(Point3d)*vs.size());
    mesh_ = msh;

    PRINT_MSG("Update normals\n");
    update_normals();
}

void TetViewerFrame::modes_params()
{
    if ( !mesh_ || maxMAmp_.empty() || modes_[0].size() != mesh_->num_vertices()*3 ) return;

    modesVisDiag_.show();
    canvas->set_disp_mode(TetViewerCanvas::SHOW_MODES);
}

void TetViewerFrame::set_mode_sel(int mId)
{
    modeSel_ = mId;
    canvas->update();
}

void TetViewerFrame::mode_scale_changed(double s)
{
    modeScale_ = s;
    canvas->update();
}

/*
 * Load vibration modes
 */
void TetViewerFrame::load_modes()
{
    if ( !mesh_ ) return;
    QString file = QFileDialog::getOpenFileName(this,
            "Select the mesh file", ".",
            "modal file (*.ev);;All file (*)");
    if ( file.isEmpty() ) return;

#if QT_VERSION >= 0x050000
    ifstream fin(file.toLatin1().data(), ios::binary);
    if ( fin.fail() )
    {
        cerr << "Fail to open the file: " << file.toLatin1().data() << endl;
        return;
    }
#else
    ifstream fin(file.toAscii().data(), ios::binary);
    if ( fin.fail() )
    {
        cerr << "Fail to open the file: " << file.toAscii().data() << endl;
        return;
    }
#endif

    int n, m;
    fin.read((char *)&m, sizeof(int));  // size of eigen-problem
    fin.read((char *)&n, sizeof(int));  // # of eigen value

    cerr << "INFO: load " << n << "  " << m << endl;
    if ( m != (int)mesh_->num_vertices()*3 )
    {
        cerr << "ERROR: the size of eigenvectors cannot match the number of vertices" << endl;
        return;
    }

    modes_.resize(n);
    maxMAmp_.resize(n);
    fin.read((char*)maxMAmp_.data(), sizeof(double)*n);
    for(int i = 0;i < n;++ i)
    {
        modes_[i].resize(m);
        fin.read((char*)modes_[i].data(), sizeof(double)*m);

        double mv = 0;
        for(int j = 0;j < m;j += 3)
            mv = fmax(mv, sqrt(M_SQR(modes_[i][j])+M_SQR(modes_[i][j+1])+M_SQR(modes_[i][j+2])));
        maxMAmp_[i] = mv;
    }
    if ( fin.fail() )
    {
#if QT_VERSION >= 0x050000
        cerr << "ERROR: Fail to read modes file: " << file.toLatin1().data() << endl;
#else
        cerr << "ERROR: Fail to read modes file: " << file.toAscii().data() << endl;
#endif
        return;
    }
    fin.close();

    cout << "Load " << n << " modes successfully" << endl;
    assert(n > 0);
    modesVisDiag_.set_mode_num(n);
}

void TetViewerFrame::update_normals()
{
    nml_.resize(mesh_->num_vertices());
    const vector<Tuple3ui>& tgl = mesh_->surface_indices();
    vector<Vector3d> tglnmls(tgl.size());
    for(size_t i = 0;i < tgl.size();++ i)
    {
        tglnmls[i] = Triangle<double>::normal(
                vtx_[tgl[i][0]], vtx_[tgl[i][1]], vtx_[tgl[i][2]]);
        if ( tglnmls[i].length_sqr() < 1E-24 )
        {
            PRINT_ERROR("triangle has zero area: %.30lf\n",
                    tglnmls[i].length_sqr());
            exit(1);
        }
        tglnmls[i].normalize();
    }

    memset(&nml_[0], 0, sizeof(double)*nml_.size());
    for(size_t i = 0;i < tgl.size();++ i)
    {
        const Vector3d& n = tglnmls[i];
        nml_[tgl[i][0]] += n * Triangle<double>::angle(
                vtx_[tgl[i][2]], vtx_[tgl[i][0]], vtx_[tgl[i][1]]);
        nml_[tgl[i][1]] += n * Triangle<double>::angle(
                vtx_[tgl[i][0]], vtx_[tgl[i][1]], vtx_[tgl[i][2]]);
        nml_[tgl[i][2]] += n * Triangle<double>::angle(
                vtx_[tgl[i][1]], vtx_[tgl[i][2]], vtx_[tgl[i][0]]);
    }

    for(size_t i = 0;i < nml_.size();++ i) 
        if ( nml_[i].length_sqr() > 1E-14 ) nml_[i].normalize();
}

void TetViewerFrame::check_useless_vtx()
{
    if ( !mesh_ ) 
    {
        PRINT_WARNING("No mesh is loaded yet\n");
        return;
    }

    bool * used = new bool[ mesh_->num_vertices() ];
    memset(used, false, sizeof(bool)*mesh_->num_vertices());

    const vector<TetMesh<double>::TetIdx>& idx = mesh_->tet_indices();
    for(size_t i = 0;i < idx.size();++ i)
    {
        used[idx[i][0]] = true;
        used[idx[i][1]] = true;
        used[idx[i][2]] = true;
        used[idx[i][3]] = true;
    }

    for(size_t i = 0;i < mesh_->num_vertices();++ i)
        if ( !used[i] ) 
        {
            PRINT_WARNING("Free vertex that is not in use is detected! VID=%d\n", (int)i);
            return;
        }
    PRINT_MSG("No Free Vertex detected\n");
}

// ============================================================================
int main(int argc, char* argv[])
{
    QApplication application(argc, argv);

    TetViewerFrame mainwnd;
    mainwnd.show();

    return application.exec();
}
