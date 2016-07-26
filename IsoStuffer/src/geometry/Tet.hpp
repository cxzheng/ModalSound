/******************************************************************************
 *  File: Tet.hpp
 *  Single tetrahedron in 3D space
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
#ifndef GEOMETRY_TETRAHEDRON_HPP
#   define GEOMETRY_TETRAHEDRON_HPP

#include <stdlib.h>
#include <stdio.h>
#include "Point3.hpp"
#include "linearalgebra/Vector3.hpp"
#include "linearalgebra/Matrix3.hpp"
#include "geometry/Triangle.hpp"
#include "utils/math.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

template <typename T>
class Tet
{
    public:
        typedef Point3<T>   TVtx;

        /* =========== constructors ========== */
        Tet(TVtx* v0, TVtx* v1, TVtx* v2, TVtx* v3)
        {
            m_vertices[0] = v0;
            m_vertices[1] = v1;
            m_vertices[2] = v2;
            m_vertices[3] = v3;

            init();
        }
        Tet(TVtx* v[])
        {
            m_vertices[0] = v[0];
            m_vertices[1] = v[1];
            m_vertices[2] = v[2];
            m_vertices[3] = v[3];

            init();
        }

        Tet() { }

        void init(TVtx* v0, TVtx* v1, TVtx* v2, TVtx* v3)
        {
            m_vertices[0] = v0;
            m_vertices[1] = v1;
            m_vertices[2] = v2;
            m_vertices[3] = v3;

            init();
        }

        /*
         * update the vertex pointer.
         * If the vertex array of tet mesh is changed, this method should be called
         * to make sure all the old tet in that tet mesh maintains the correct 
         * vertex points.
         */
        inline void conform_vtx_pointers(TVtx* v0, TVtx* v1, TVtx* v2, TVtx* v3)
        {
            m_vertices[0] = v0;
            m_vertices[1] = v1;
            m_vertices[2] = v2;
            m_vertices[3] = v3;

            m_faces[0].init(m_vertices[0], m_vertices[1], m_vertices[3]);
            m_faces[1].init(m_vertices[0], m_vertices[2], m_vertices[1]);
            m_faces[2].init(m_vertices[0], m_vertices[3], m_vertices[2]);
            m_faces[3].init(m_vertices[1], m_vertices[2], m_vertices[3]);
        }

        /*
         * Solid angle subtended by i-th vertex
         */
        T solid_angle(int iv) const;

        /* =========== Retrival methods ============ */
        TVtx** vertices() { return m_vertices; }

        TVtx* vtx(int n)
        {
            assert(n >= 0 && n < 4);
            return m_vertices[n];
        }

        const TVtx* vtx(int n) const
        {
            assert(n >= 0 && n < 4);
            return m_vertices[n];
        }

        const TVtx& operator [] (int n) const
        {
            assert(n >= 0 && n < 4);
            return *m_vertices[n];
        }

        const Matrix3<T>& inverse_Dm() const
        {   return m_invDm; }

        void deformation_gradient(Matrix3<T>& mat) const;
        const Matrix3<T>& deformation_gradient();       //  $$TESTED
        /*! return weighted normals in material space */
        const Vector3<T>* weighted_normals() const 
        { return m_b; }

        /*! return current volume of this tetrahedron */
        T volume() const;           //                      $$TESTED
        T signed_volume() const;    // right-hand system    $$TESTED

        static inline T volume(const Point3<T>& p1, const Point3<T>& p2,
                const Point3<T>& p3, const Point3<T>& p4)
        {
            // formula for a tet volume with vertices (a,b,c,d) is:
            // |(a - d) dot ((b - d) cross (c - d))| / 6
            const Vector3<T> ad = p2 - p1;
            const Vector3<T> bd = p3 - p1;
            const Vector3<T> cd = p4 - p1;

            return M_ABS(ad.dotProduct(bd.crossProduct(cd))) / (T)6;
        }

        /*! get the geometric center */
        Point3<T> geometric_center() const
        {
            return (*m_vertices[0] + *m_vertices[1] + *m_vertices[2] + *m_vertices[3])*0.25;
        }

    private:
        void init();

    private:
        TVtx*       m_vertices[4];
        Vector3<T>  m_b[4];     // weighted normals in material space
        Triangle<T> m_faces[4];
        Matrix3<T>  m_invDm;    // inverse of D_m
        Matrix3<T>  m_F;        // deformation gradient
};

/*
 * the vertex layout is assumed to be like this:
 *
 *            v3
 *           /| \
 *          / |  \
 *         /  |   \  
 *        /f0 |    \ 
 *       /    |     \
 *      /     |      \  
 *     /      |       \
 *    v0- - - - - - - -v2
 *     \      |       /
 *       \    | f1  /  
 *         \  |   / 
 *           \| / 
 *            v1
 */  
template <typename T>
void Tet<T>::init() 
{
    m_faces[0] = Triangle<T>(m_vertices[0], m_vertices[1], m_vertices[3]);
    m_faces[1] = Triangle<T>(m_vertices[0], m_vertices[2], m_vertices[1]);
    m_faces[2] = Triangle<T>(m_vertices[0], m_vertices[3], m_vertices[2]);
    m_faces[3] = Triangle<T>(m_vertices[1], m_vertices[2], m_vertices[3]);

    // compute weighted normals
    Vector3<T> wn[4];
    wn[0] = m_faces[0].weighted_normal();
    wn[1] = m_faces[1].weighted_normal();
    wn[2] = m_faces[2].weighted_normal();
    wn[3] = m_faces[3].weighted_normal();

    const T ss = (T)(-1.) / (T)3;
    m_b[0] = (wn[0] + wn[1] + wn[2]) * ss;
    m_b[1] = (wn[0] + wn[1] + wn[3]) * ss;
    m_b[2] = (wn[3] + wn[1] + wn[2]) * ss;
    m_b[3] = (wn[0] + wn[3] + wn[2]) * ss;

    // compute inverse of Dm
    Matrix3<T>((*m_vertices[1]) - (*m_vertices[0]),
               (*m_vertices[2]) - (*m_vertices[0]),
               (*m_vertices[3]) - (*m_vertices[0])).inverse(m_invDm);

    if ( signed_volume() < 1E-16 )
    {
        fprintf(stderr, "ERROR: Tetrahedron configuration error! V = %.26lf\n",
                signed_volume());
        exit(1);
    }
}

template <typename T>
T Tet<T>::signed_volume() const
{
    Vector3<T> ad = (*m_vertices[1]) - (*m_vertices[0]);    // col[0]
    Vector3<T> bd = (*m_vertices[2]) - (*m_vertices[0]);    // col[1]
    Vector3<T> cd = (*m_vertices[3]) - (*m_vertices[0]);    // col[2]

    // (col[1] x col[2]) . col[0]
    return ad.dotProduct(bd.crossProduct(cd)) / (T)6;
}

template <typename T>
T Tet<T>::volume() const
{
    // formula for a tet volume with vertices (a,b,c,d) is:
    // |(a - d) dot ((b - d) cross (c - d))| / 6
    Vector3<T> ad = (*m_vertices[1]) - (*m_vertices[0]);
    Vector3<T> bd = (*m_vertices[2]) - (*m_vertices[0]);
    Vector3<T> cd = (*m_vertices[3]) - (*m_vertices[0]);

    return M_ABS(ad.dotProduct(bd.crossProduct(cd))) / (T)6;
}

template <typename T>
const Matrix3<T>& Tet<T>::deformation_gradient()
{
    m_F = Matrix3<T>(*m_vertices[1] - *m_vertices[0],
                     *m_vertices[2] - *m_vertices[0],
                     *m_vertices[3] - *m_vertices[0]) * m_invDm;

    return m_F;
}

template <typename T>
void Tet<T>::deformation_gradient(Matrix3<T>& mat) const
{
    mat = Matrix3<T>(*m_vertices[1] - *m_vertices[0],
                     *m_vertices[2] - *m_vertices[0],
                     *m_vertices[3] - *m_vertices[0]) * m_invDm;
}

/*
 * compute the solid angle subtended by the i-th vertex 
 *
 * The formula to computing the solid angle is given in the paper
 * A. Van Oosterom et. al. The Solid Angle of a Plane Triangle, IEEE
 * Transactions on Biomedical Engineering, Vol. BME-30, No. 2, Feb. 1983
 *
 * NOTE: The returned value should always be positive
 */
template <typename T>
T Tet<T>::solid_angle(int iv) const
{
    assert(iv >= 0 && iv < 4);
    static const int IDX[][3] = { {2, 3, 1}, {0, 3, 2}, 
                                  {1, 3, 0}, {2, 1, 0} };
    Matrix3<T> mat(*m_vertices[IDX[iv][0]] - *m_vertices[iv],
                   *m_vertices[IDX[iv][1]] - *m_vertices[iv],
                   *m_vertices[IDX[iv][2]] - *m_vertices[iv]);

    T vlen0 = mat.cols[0].length();
    T vlen1 = mat.cols[1].length();
    T vlen2 = mat.cols[2].length();

    return 2.*atan2(mat.det(), 
            vlen0*vlen1*vlen2 + 
            mat.cols[0].dotProduct(mat.cols[1]) * vlen2 +
            mat.cols[0].dotProduct(mat.cols[2]) * vlen1 + 
            mat.cols[1].dotProduct(mat.cols[2]) * vlen0);
}

#ifdef USE_NAMESPACE
}
#endif

#endif
