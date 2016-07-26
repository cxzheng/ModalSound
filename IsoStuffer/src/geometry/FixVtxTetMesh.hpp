/******************************************************************************
 *  File: FixVtxTetMesh.hpp
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
#ifndef GEOMETRY_FIXVTXTETMESH_HPP
#   define GEOMETRY_FIXVTXTETMESH_HPP

#include "geometry/TetMesh.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

/*
 * The tet mesh with fixed vertices. It assumes the FIRST N vertices
 * in the m_vertices list are fixed.
 */
template <typename T>
class FixVtxTetMesh : public TetMesh<T>
{
    using TetMesh<T>::m_vertices;
    using TetMesh<T>::m_restPos;

    public:
        FixVtxTetMesh():m_numFixed(0) { }

        FixVtxTetMesh(const FixVtxTetMesh<T>& rhs):TetMesh<T>(rhs),
                m_numFixed(rhs.m_numFixed)
        {}

        /*!
         * set the number of fixed vertices, N
         * Assume the FIRST N vertices in m_vertices are fixed
         */
        void set_fixed_vtx(int n)
        {
            m_numFixed = n;
            assert(m_numFixed <= m_vertices.size());
        }

        /*!
         * Add vertex into the current tet mesh.
         * Note: if the vertex is fixed, then the vertex is move to the beginning
         *
         * NOTE that the method can only be called before any tet is added into this
         *      tet mesh.
         *
         * NOTE the returned value is the vertex ID of the added vertex
         */
        template <typename FromT>
        int add_fixed_vertex(const Point3<FromT>& rest, const Point3<FromT>& vnow, bool fixed)
        {
            if ( fixed )
            {
                if ( m_numFixed == (int)m_vertices.size() )  // all the existing vertices are fixed
                {
                    m_restPos.push_back(rest);
                    m_vertices.push_back(vnow);
                }
                else
                {
                    m_restPos.push_back(m_restPos[m_numFixed]);
                    m_vertices.push_back(m_vertices[m_numFixed]);

                    m_restPos[m_numFixed] = rest;
                    m_vertices[m_numFixed] = vnow;
                }
                ++ m_numFixed;
                return m_numFixed - 1;
            }
            else
            {
                m_restPos.push_back(rest);
                m_vertices.push_back(vnow);

                return m_vertices.size() - 1;
            }
        }

        bool is_fixed_vertex(int vid) const
        {  
            assert(vid >= 0);
            return vid < m_numFixed; 
        }

        size_t num_fixed_vertices() const 
        {  return m_numFixed; }

        size_t num_free_vertices() const 
        {  return m_vertices.size() - (size_t)m_numFixed; }

    private:
        int     m_numFixed;
};

#ifdef USE_NAMESPACE
}
#endif
#endif

