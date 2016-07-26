/******************************************************************************
 *  File: DistValProc.cpp
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
#include "DistValProc.h"
#include <string.h>

DistValProc::DistValProc(const TriangleMesh<double>* mesh, 
        int res, int nlevel):mp_mesh(mesh), m_obbtree(mesh)
{
    // find the min/max point of mesh
    m_minP.set(1E+100, 1E+100, 1E+100);
    m_maxP.set(-1E+100, -1E+100, -1E+100);
    const std::vector<Point3d>& vtx = mp_mesh->vertices();
    for(size_t i = 0;i < vtx.size();++ i)
    {
        m_minP.x = fmin(m_minP.x, vtx[i].x);
        m_minP.y = fmin(m_minP.y, vtx[i].y);
        m_minP.z = fmin(m_minP.z, vtx[i].z);

        m_maxP.x = fmax(m_maxP.x, vtx[i].x);
        m_maxP.y = fmax(m_maxP.y, vtx[i].y);
        m_maxP.z = fmax(m_maxP.z, vtx[i].z);
    }

    update_mesh_params(res, nlevel);

    //// pre-compute pseudo-normals
    compute_tgl_pseudo_normals();
    compute_edge_pseudo_normals();
    compute_vtx_pseudo_normals();
    compute_tgl_collision_info();
}

void DistValProc::update_mesh_params(int res, int nlevel, int margin)
{
    if ( res < (1<<(nlevel-1)) )
    {
        std::cerr << "ERROR: nLevel is too large for the resolution" << std::endl;
        exit(1);
    }

    m_minCorner = m_maxP - m_minP;      // scale of the object
    double L = fmax(m_minCorner.z, fmax(m_minCorner.x, m_minCorner.y));
    const int NL = 1 << nlevel;

    m_gds = L / (double)(res-margin);
    m_res.set(res, res, res);
    for(int i = 0;i < 3;++ i)
    {
        while ( (m_res[i] / 2)*m_gds > m_minCorner[i] + margin*m_gds &&
                m_res[i] >= NL ) m_res[i] = m_res[i] / 2;
    }
    Tuple3d llen(m_res);
    llen *= m_gds;
    m_minCorner = (m_minP + m_maxP - llen)*0.5;

    std::cout << "INFO: resolution: " << m_res << "  # of level: " << nlevel << std::endl;
    std::cout << "      max point: " << m_maxP << std::endl;
    std::cout << "      min point: " << m_minP << std::endl;
    std::cout << "      Grid Size: " << m_gds << std::endl;
    std::cout << "      Min point: " << m_minCorner << std::endl;
}

void DistValProc::compute_tgl_collision_info()
{
    const std::vector<Point3d>&  vtx = mp_mesh->vertices(); 
    const std::vector<Tuple3ui>& tgl = mp_mesh->surface_indices();

    m_tglCollInfo.resize(tgl.size());
    for(size_t i = 0;i < tgl.size();++ i)
    for(int j = 0;j < 3;++ j)
    {
        int vid0 = tgl[i][j], vid1 = tgl[i][(j+1)%3], vid2 = tgl[i][(j+2)%3];
        m_tglCollInfo[i].vtxId[j] = vid0;
        m_tglCollInfo[i].sidelen[j] = vtx[vid0].distance(vtx[vid1]);

        Vector3d vec01 = vtx[vid1] - vtx[vid0];
        m_tglCollInfo[i].QT[j].cols[0] = vec01 / m_tglCollInfo[i].sidelen[j];
        m_tglCollInfo[i].QT[j].cols[2] = vec01.crossProduct(vtx[vid2] - vtx[vid0]);
        m_tglCollInfo[i].QT[j].cols[2].normalize();
        m_tglCollInfo[i].QT[j].cols[1] = m_tglCollInfo[i].QT[j].cols[2].crossProduct(
                m_tglCollInfo[i].QT[j].cols[0]);

        m_tglCollInfo[i].Q[j] = m_tglCollInfo[i].QT[j].transpose();
        m_tglCollInfo[i].x0[j] = -(m_tglCollInfo[i].Q[j] * vtx[vid0]);
    }
}

void DistValProc::compute_tgl_pseudo_normals()
{
    const std::vector<Point3d>&   vtx = mp_mesh->vertices(); 
    const std::vector<Tuple3ui>&  tgl = mp_mesh->surface_indices();

    m_tglPseudoNml.resize(tgl.size());
    for(size_t i = 0;i < tgl.size();++ i)
    {
        m_tglPseudoNml[i] = Triangle<double>::normal(
                vtx[tgl[i][0]], vtx[tgl[i][1]], vtx[tgl[i][2]]);
        if ( m_tglPseudoNml[i].lengthSqr() < 1E-24 )
        {
            fprintf(stderr, "ERROR: triangle has zero area: (%.10f %.10f %.10f)(%.10f %.10f %.10f)(%.10f %.10f %.10f) --> %.30lf\n",
                    vtx[tgl[i][0]].x, vtx[tgl[i][0]].y, vtx[tgl[i][0]].z,
                    vtx[tgl[i][1]].x, vtx[tgl[i][1]].y, vtx[tgl[i][1]].z,
                    vtx[tgl[i][2]].x, vtx[tgl[i][2]].y, vtx[tgl[i][2]].z,
                    m_tglPseudoNml[i].lengthSqr());
            exit(1);
        }

        /////// Bug fix according to Tim L
        m_tglPseudoNml[i].normalize();
    }
}

void DistValProc::compute_edge_pseudo_normals()
{
    m_edgePseudoNml.clear();

    const std::vector<Tuple3ui>&  tgl = mp_mesh->surface_indices();
    for(size_t i = 0;i < tgl.size();++ i)
    {
        int a, b;
        // edge 0
        a = std::min(tgl[i][0], tgl[i][1]);
        b = std::max(tgl[i][0], tgl[i][1]);
        m_edgePseudoNml[a][b] += m_tglPseudoNml[i];

        // edge 1
        a = std::min(tgl[i][1], tgl[i][2]);
        b = std::max(tgl[i][1], tgl[i][2]);
        m_edgePseudoNml[a][b] += m_tglPseudoNml[i];

        // edge 1
        a = std::min(tgl[i][2], tgl[i][0]);
        b = std::max(tgl[i][2], tgl[i][0]);
        m_edgePseudoNml[a][b] += m_tglPseudoNml[i];
    }
}

void DistValProc::compute_vtx_pseudo_normals()
{
    const std::vector<Point3d>&  vtx = mp_mesh->vertices(); 
    const std::vector<Tuple3ui>& tgl = mp_mesh->surface_indices();

    m_vtxPseudoNml.resize(vtx.size());
    memset(&m_vtxPseudoNml[0], 0, sizeof(Vector3d)*m_vtxPseudoNml.size());

    for(size_t i = 0;i < tgl.size();++ i)
    {
        const Vector3d& nml = m_tglPseudoNml[i];

        m_vtxPseudoNml[tgl[i][0]] += nml * Triangle<double>::angle(
                vtx[tgl[i][2]], vtx[tgl[i][0]], vtx[tgl[i][1]]);
        m_vtxPseudoNml[tgl[i][1]] += nml * Triangle<double>::angle(
                vtx[tgl[i][0]], vtx[tgl[i][1]], vtx[tgl[i][2]]);
        m_vtxPseudoNml[tgl[i][2]] += nml * Triangle<double>::angle(
                vtx[tgl[i][1]], vtx[tgl[i][2]], vtx[tgl[i][0]]);
    }
}
