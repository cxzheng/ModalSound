#ifndef GEOMETRY_FIXEDVTXTETMESH_HPP
#   define GEOMETRY_FIXEDVTXTETMESH_HPP

#include "geometry/TetMesh.hpp"

#ifdef USE_NAMESPACE
namespace sploosh
{
#endif

/*
 * The tet mesh with fixed vertices. It assumes the FIRST N vertices
 * in the vertices_ list are fixed.
 */
template <typename T>
class FixedVtxTetMesh : public TetMesh<T>
{
    using TetMesh<T>::vertices_;
    using TetMesh<T>::restPos_;

    public:
        FixedVtxTetMesh():numFixed_(0) { }

        FixedVtxTetMesh(const FixedVtxTetMesh<T>& rhs):TetMesh<T>(rhs),
                numFixed_(rhs.numFixed_)
        {}

        /*!
         * set the number of fixed vertices, N
         * Assume the FIRST N vertices in vertices_ are fixed
         */
        void set_fixed_vtx(int n)
        {
            numFixed_ = n;
            assert(numFixed_ <= vertices_.size());
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
                if ( numFixed_ == (int)vertices_.size() )  // all the existing vertices are fixed
                {
                    restPos_.push_back(rest);
                    vertices_.push_back(vnow);
                }
                else
                {
                    restPos_.push_back(restPos_[numFixed_]);
                    vertices_.push_back(vertices_[numFixed_]);

                    restPos_[numFixed_]  = rest;
                    vertices_[numFixed_] = vnow;
                }
                ++ numFixed_;
                return numFixed_ - 1;
            }
            else
            {
                restPos_.push_back(rest);
                vertices_.push_back(vnow);

                return vertices_.size() - 1;
            }
        }

        bool is_fixed_vertex(int vid) const
        {  
            assert(vid >= 0);
            return vid < numFixed_; 
        }

        size_t num_fixed_vertices() const 
        {  return numFixed_; }

        size_t num_free_vertices() const 
        {  return vertices_.size() - (size_t)numFixed_; }

    private:
        int     numFixed_;
};

#ifdef USE_NAMESPACE
}
#endif
#endif

