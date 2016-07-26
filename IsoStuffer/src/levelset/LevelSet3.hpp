/******************************************************************************
 *  File: LevelSet3.hpp
 *  A discrete level-set domain divided by voxels in 3D space
 *  Copyright (c) 2007 by Changxi Zheng
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
#ifndef LEVEL_SET_3D_HPP
#   define LEVEL_SET_3D_HPP

#include "constants.h"
#include "utils/arrays.hpp"
#include "utils/trilinear.hpp"
#include "linearalgebra/Vector3.hpp"
#include "geometry/Point3.hpp"

template <typename T, class DistFunc, bool BeLazy>
class LevelSet3
{
    public:
        LevelSet3(const Point3<T>& minPt, const Point3<T>& maxPt, DistFunc func);

        /*
         * given an point in the space, return the signed distance value
         */
        T distance(const Point3<T>& pt);

        /*
         * Given a position, check the distance value at that point.
         * if the distance is negative (inside of the object), return 
         * that distance and the normal direction at that point.
         */
        T negative_dist_with_normal(const Point3<T>& pt, Vector3<T>& nml);

        /*
         * get the distance field value at the voxel (cx,cy,cz)
         */
        T levelset_value(int cx, int cy, int cz);

        void get_normal(int cx, int cy, int cz, Vector3<T>& nml);

        void voxel_coord(const Point3<T>& pos, Point3i& coord, Tuple3<T>& alpha) const;

        bool inside_box(const Point3<T>& pt) const
        {
            return (pt.x > m_minBound.x && pt.x < m_maxBound.x &&
                    pt.y > m_minBound.y && pt.y < m_maxBound.y &&
                    pt.z > m_minBound.z && pt.z < m_maxBound.z);
        }

    private:
        typedef boost::multi_array<T, 3>            ScalarField;
        typedef boost::multi_array<Vector3<T>, 3>   VectorField;

        ScalarField     m_distfield;    // discrete distance field
        VectorField     m_nmlfield;     // discrete normal field

        Point3<T>       m_minPt;        // the minimum corner, described by this level-set
        Point3<T>       m_maxPt;        // the maximum corner, described by this level-set
        Point3<T>       m_minBound;
        Point3<T>       m_maxBound;

        /* 
         * size of each voxel
         * each voxel has the same size on three dimensions
         */
        REAL            m_VS;
        REAL            m_invVS;
        /*
         * The resolution of this level-set grids NX x NY x NZ
         * there are (NX+1) points along X-axis, so do Y- and Z- axis.
         */
        Tuple3ui        m_res;          // resolutions
        DistFunc        _dist_func;     // function to compute the distance value
};

///////////////////////////////////////////////////////////////////////////////

template <typename T, class DistFunc, bool BeLazy>
LevelSet3<T, DistFunc, BeLazy>::LevelSet3(
        const Point3<T>& minPt, const Point3<T>& maxPt, DistFunc func):
        m_minPt(minPt), m_maxPt(maxPt), _dist_func(func)
{
    printf("MSG: Allocate LevelSet Voxels ... ");
    Vector3<T> dd = m_maxPt - m_minPt;
    m_VS = fmin(dd.z*0.5, fmin(dd.y*0.5, fmin(LEVELSET_VOXEL_SIZE, dd.x*0.5)));
    m_invVS = 1. / m_VS;

    m_res.set((int)(dd.x * m_invVS) + 9,
              (int)(dd.y * m_invVS) + 9,
              (int)(dd.z * m_invVS) + 9);

    Vector3<T> dd2(m_res.x * m_VS, m_res.y * m_VS, m_res.z * m_VS);
    dd2 -= dd;
    dd2 *= 0.5;
    m_minPt -= dd2;
    m_maxPt += dd2;

    m_minBound = m_minPt + 0.5*m_VS;
    m_maxBound = m_maxPt - 0.5*m_VS;
    //// init boost multi-array
    m_distfield.resize(boost::extents[m_res.z+1][m_res.y+1][m_res.x+1]);
    m_nmlfield.resize(boost::extents[m_res.z+1][m_res.y+1][m_res.x+1]);

    zero_multi_array(m_distfield);
    zero_multi_array(m_nmlfield);
    printf("[OK]\n");
}

template <typename T, class DistFunc, bool BeLazy>
void LevelSet3<T, DistFunc, BeLazy>::voxel_coord(
        const Point3<T>& pos, Point3i& coord, Tuple3<T>& alpha) const
{
    alpha.set((pos.x - m_minPt.x)*m_invVS, 
              (pos.y - m_minPt.y)*m_invVS, 
              (pos.z - m_minPt.z)*m_invVS);
    coord.set((int)alpha.x, (int)alpha.y, (int)alpha.z);
    alpha -= coord;
}

template <typename T, class DistFunc, bool BeLazy>
T LevelSet3<T, DistFunc, BeLazy>::distance(const Point3<T>& pt)
{

    if ( !inside_box(pt) ) return 1E+10;

    Point3i      coord;
    Tuple3<REAL> alpha;
    voxel_coord(pt, coord, alpha);
    T iso[2][2][2];

    for(int iz = 0;iz < 2;++ iz)
    for(int iy = 0;iy < 2;++ iy)
    for(int ix = 0;ix < 2;++ ix)
    {
        const int cx = coord.x + ix;
        const int cy = coord.y + iy;
        const int cz = coord.z + iz;
        iso[iz][iy][ix] = levelset_value(cx, cy, cz);
    }
    T ret = trilinear_interpolate((const T (*)[2][2])iso, (const T*)(&alpha));
    return ret;
}

/*
 * - get the voxel coordinate of the given point
 * - get the level-set value at 8 point around that pt
 *
 * NOTE: pt should be given in the INITIAL configuration of the rigid body
 *       because LevelSet box is created based on the initial state of the rigid
 *       body.
 *       the output normal direction is therefore also in the rigid body's
 *       INITIAL configuration.
 *
 * \param de is the argument that is passed to _dist_func when computing
 *           the distance value at the given point if necessary. Because we are
 *           using OBB to find nearest triangles to the given point, the given
 *           argument is the estimation of the distance from the point to the 
 *           object associated with this level set.
 *        de is the distance estimation, no square here.
 */
template <typename T, class DistFunc, bool BeLazy>
T LevelSet3<T, DistFunc, BeLazy>::negative_dist_with_normal(
        const Point3<T>& pt, Vector3<T>& nml)
{
    if ( !inside_box(pt) ) return 1E+10;

    Point3i      coord;
    Tuple3<REAL> alpha;
    voxel_coord(pt, coord, alpha);
    T iso[2][2][2];

    for(int iz = 0;iz < 2;++ iz)
    for(int iy = 0;iy < 2;++ iy)
    for(int ix = 0;ix < 2;++ ix)
    {
        const int cx = coord.x + ix;
        const int cy = coord.y + iy;
        const int cz = coord.z + iz;
        iso[iz][iy][ix] = levelset_value(cx, cy, cz);
    }
    T ret = trilinear_interpolate((const T (*)[2][2])iso, (const T*)(&alpha));
    if ( ret > 0. ) return ret;

    //// evaluate the normal 
    Vector3<T> vecs[2][2][2];

    for(int iz = 0;iz < 2;++ iz)
    for(int iy = 0;iy < 2;++ iy)
    for(int ix = 0;ix < 2;++ ix)
    {
        const int cx = coord.x + ix;
        const int cy = coord.y + iy;
        const int cz = coord.z + iz;
        get_normal(cx, cy, cz, vecs[iz][iy][ix]);
    }

    nml = trilinear_interpolate((const Vector3<T> (*)[2][2])vecs, (const T*)&alpha);
    nml.normalize();
    return ret;
}

/*
 * \param dse Estimate of the square of the distance
 */
template <typename T, class DistFunc, bool BeLazy>
T LevelSet3<T, DistFunc, BeLazy>::levelset_value(int cx, int cy, int cz)
{
    if ( BeLazy && m_distfield[cz][cy][cx] == 0.0 ) // assuming no point has exactly 0 distance value
        return (m_distfield[cz][cy][cx] = _dist_func(Point3<REAL>(
                        cx*m_VS + m_minPt.x, cy*m_VS + m_minPt.y, 
                        cz*m_VS + m_minPt.z)));
    else
        return m_distfield[cz][cy][cx];
}

/*
 * \param de Estimate of the distance
 */
template <typename T, class DistFunc, bool BeLazy>
void LevelSet3<T, DistFunc, BeLazy>::get_normal(
        int cx, int cy, int cz, Vector3<T>& nml)
{
    if ( !BeLazy || m_nmlfield[cz][cy][cx].lengthSqr() != 0. )
    {
        nml = m_nmlfield[cz][cy][cx];
        return;
    }

    //// interpolate the normal at (cx,cy,cz)
    const int xplus  = std::min(cx+1, (int)m_res.x);
    const int xminus = std::max(0, cx-1);
    const int yplus  = std::min(cy+1, (int)m_res.y);
    const int yminus = std::max(0, cy-1);
    const int zplus  = std::min(cz+1, (int)m_res.z);
    const int zminus = std::max(0, cz-1);

    m_nmlfield[cz][cy][cx].set(
            (levelset_value(xplus, cy, cz ) - 
             levelset_value(xminus, cy, cz)) / (m_VS*(xplus - xminus)),
            (levelset_value(cx, yplus, cz ) - 
             levelset_value(cx, yminus, cz)) / (m_VS*(yplus - yminus)),
            (levelset_value(cx, cy, zplus ) - 
             levelset_value(cx, cy, zminus)) / (m_VS*(zplus - zminus)));
    m_nmlfield[cz][cy][cx].normalize();

    nml = m_nmlfield[cz][cy][cx];
}

#endif
