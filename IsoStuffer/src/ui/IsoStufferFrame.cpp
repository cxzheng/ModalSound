/******************************************************************************
 *  File: IsoStufferFrame.cpp
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
#include "IsoStufferFrame.h"
#include <QFileDialog>
#include <QCloseEvent>

#include <algorithm>
#include <map>
#include "io/TglMeshReader.hpp"
#include "io/TetMeshWriter.hpp"
#include "utils/macros.h"

void IsoStufferFrame::closeEvent(QCloseEvent *event)
{
    params_.close();
    event->accept();
}

void IsoStufferFrame::convert_to_tetmesh()
{
    if ( !dvproc_ ) return;

    delete pstuffer_;
    dvproc_->update_mesh_params(params_.resolution(), 
            params_.num_levels(), params_.margin());
    pstuffer_ = new TIsoStuffer(dvproc_->resolution(),
            params_.num_levels(), dvproc_->grid_size(),
            dvproc_->min_point(), params_.alpha_long(),
            params_.alpha_short(), *pdistfunc_);

    delete pmesh_;
    pmesh_ = new TetMesh<double>;

    pstuffer_->create_mesh(pmesh_);
    pmesh_->update_surface();
    curStep_ = 6;
    canvas->update();
}

void IsoStufferFrame::next_step()
{
    switch ( curStep_ )
    {
        case 0: // create boundary octants
            if ( !dvproc_ ) return;

            delete pstuffer_;
            dvproc_->update_mesh_params(params_.resolution(), 
                    params_.num_levels(), params_.margin());
            canvas->scale_ = dvproc_->grid_size() * 0.15;
            pstuffer_ = new TIsoStuffer(dvproc_->resolution(),
                    params_.num_levels(), dvproc_->grid_size(),
                    dvproc_->min_point(), params_.alpha_long(),
                    params_.alpha_short(), *pdistfunc_);

            pstuffer_->create_boundary_octants();
            ++ curStep_;
            canvas->update();
            return;
        case 1: // create octree
            pstuffer_->create_octree();
            std::cerr << "CHECK CHILDREN MASK ... I" << std::endl;
            pstuffer_->oct_tree().check_children_mask();

            pstuffer_->fill_highest_level();
            //std::cerr << "CHECK CHILDREN MASK ... II" << std::endl;
            //pstuffer_->oct_tree().check_children_mask();
            ++ curStep_;
            canvas->update();
            return;
        case 2: // weak balance
            pstuffer_->weak_balance();
            std::cerr << "CHECK CHILDREN MASK ... III" << std::endl;
            pstuffer_->oct_tree().check_children_mask();
            pstuffer_->update_bitmasks();
            ++ curStep_;
            canvas->update();
            return;
        case 3:
            pstuffer_->create_background_tets();

            //// create background mesh
            delete pmesh_;
            pmesh_ = new TetMesh<double>;
            pstuffer_->background_tet_mesh(pmesh_);
            pmesh_->update_surface();

            ++ curStep_;
            canvas->update();
            return;
        case 4:
            pstuffer_->compute_cutting_pts();
            pstuffer_->wrap_violated_pts();
            ++ curStep_;
            canvas->update();
            return;
        case 5:
            printf("INFO: extract final tets...\n");
            pstuffer_->extract_final_tets();

            delete pmesh_;
            pmesh_ = new TetMesh<double>;
            pstuffer_->extract_final_mesh(pmesh_);
            pmesh_->update_surface();

            ++ curStep_;
            canvas->update();
            return;
        case 6:
            return;
        default:
            fprintf(stderr, "Unknown step ID %d\n", curStep_);
            SHOULD_NEVER_HAPPEN(1);
    }
}

void IsoStufferFrame::reset()
{
    SAFE_DELETE(ptglmesh_);
    SAFE_DELETE(pstuffer_);
    SAFE_DELETE(pmesh_);
    SAFE_DELETE(dvproc_);
    SAFE_DELETE(pdistfunc_);
    curStep_ = 0;
}

void IsoStufferFrame::load_mesh()
{
    QString file = QFileDialog::getOpenFileName(this,
            "Select the triangle mesh file", "",
            "Obj file (*.obj);;All file (*)");

    if ( file.isEmpty() ) return;

    curFile_ = file;
#if QT_VERSION >= 0x050000
    load_tgl_mesh(curFile_.toLatin1().data());
#else
    load_tgl_mesh(curFile_.toAscii().data());
#endif
}

void IsoStufferFrame::load_tgl_mesh(const char* file)
{
    if ( curFile_.isEmpty() ) return;

    reset();
    ptglmesh_ = new TriangleMesh<double>;
    MeshObjReader::read(file, *ptglmesh_, 
            actionCenterize->isChecked(), 
            actionFlipTgl->isChecked());

    dvproc_ = new DistValProc(ptglmesh_, 
            params_.resolution(), 
            params_.num_levels());
    pdistfunc_ = new TDistFunc(dvproc_);
    canvas->scale_ = dvproc_->grid_size() * 0.15;
}

void IsoStufferFrame::check_mesh()
{
    using namespace std;

    if ( !ptglmesh_ )
    {
        cerr << "WARNING: No mesh loaded." << endl;
        return;
    }

    const vector<Point3d>&  vtx = ptglmesh_->vertices();
    const vector<Tuple3ui>& tgl = ptglmesh_->surface_indices();
    map< int, map<int, int> >  edgecnt;
    double minlength = 1E+300;
    for(size_t i = 0;i < tgl.size();++ i)
    for(int j = 0;j < 3;++ j)
    {
        int v1 = tgl[i][j], v2 = tgl[i][(j+1)%3];
        if ( v1 > v2 ) swap(v1, v2);
        minlength = min(minlength, vtx[v1].distance(vtx[v2]));
        ++ edgecnt[v1][v2];
    }

    const map< int, map<int, int> >::iterator end = edgecnt.end();
    for(map< int, map<int, int> >::iterator it1 = edgecnt.begin();
            it1 != end;++ it1)
    {
        const map<int, int>::iterator end2 = it1->second.end();
        for(map<int,int>::iterator it2 = it1->second.begin(); 
                it2 != end2;++ it2)
            if ( it2->second != 2 )
            {
                cerr << "Mesh check FAILED: " << it1->first << "--" << it2->first << "  : " << it2->second << endl;
                cerr << "   minimum distance: " << minlength << endl;
                return;
            }
    }
    cerr << "Mesh check succeed" << endl;
}

void IsoStufferFrame::save_mesh() 
{
    QString file = QFileDialog::getSaveFileName(this,
            "Saved tetrahedron file", ".", "Text file (*.tet);;All file (*)");
    if ( file.isEmpty() ) return;


#if QT_VERSION >= 0x050000
    printf("INFO: write to mesh: %s ...", file.toLatin1().data());
    FV_TetMeshWriter_Double::write_mesh(file.toLatin1().data(), *pmesh_);
#else
    printf("INFO: write to mesh: %s ...", file.toAscii().data());
    FV_TetMeshWriter_Double::write_mesh(file.toAscii().data(), *pmesh_);
#endif
    printf(" OK\n");
}

void IsoStufferFrame::export_to_stellar() 
{
    /*
    QString file = QFileDialog::getSaveFileName(this,
            "Saved tetrahedron file", ".", "All file (*)");
    if ( file.isEmpty() ) return;

    printf("INFO: export to stellar mesh: %s ...", file.toAscii().data());
    StellarMeshWriter::write_mesh(file.toAscii().data(), *pmesh_);
    printf(" OK\n");
    */
}
